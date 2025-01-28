# Load necessary libraries
import sys
import pandas as pd
import anndata as ad
import numpy as np
import scipy.sparse as sparse

## Should be able to store very large sizes:
## https://github.com/scverse/scanpy_usage/tree/master/170522_visualizing_one_million_cells

def create_anndata_from_dfs(expression_df, metadata_df, spatial_df=None, image=None, image_key = "image"):
    """
    Creates an AnnData object from Pandas DataFrames.

    Args:
        expression_df: DataFrame containing gene expression data (rows: cells, cols: genes).
        metadata_df: DataFrame containing cell metadata (rows: cells, cols: metadata).
        spatial_df (optional): DataFrame containing spatial coordinates (rows: cells, cols: X, Y).
        image (optional): numpy array containing the image
        image_key (optional): key to save the image to adata.uns or adata.obsm. Default is "image".

    Returns:
        An AnnData object or None if there is an error.
    """

    try:
        # Check for matching indices
        if not expression_df.index.equals(metadata_df.index):
            raise ValueError("Expression and metadata DataFrames must have matching indices.")
        if spatial_df is not None and not expression_df.index.equals(spatial_df.index):
            raise ValueError("Expression and spatial DataFrames must have matching indices.")

        # Create AnnData object
        adata = ad.AnnData(expression_df)

        # Add metadata to .obs
        adata.obs = metadata_df

        # Add spatial coordinates to .obsm
        if spatial_df is not None:
            adata.obsm["spatial"] = spatial_df[["X", "Y"]].values

        #Add image to .uns or .obsm
        if image is not None:
            if image.ndim == 2 or image.ndim == 3: #check if image is valid
                if adata.n_obs == image.shape[0]:
                    adata.obsm[image_key] = image
                    print("Image saved in .obsm")
                else:
                    adata.uns[image_key] = image
                    print("Image saved in .uns")
            else:
                raise ValueError("Image must be a 2D (grayscale) or 3D (RGB) numpy array")

        return adata

    except ValueError as e:
        print(f"Error creating AnnData object: {e}")
        return None

# Example usage:
# Create dummy dataframes
np.random.seed(0)
n_cells = 100
n_genes = 20
expression_data = np.random.rand(n_cells, n_genes)
expression_df = pd.DataFrame(expression_data, columns=[f"gene_{i}" for i in range(n_genes)])
expression_df.index = [f"cell_{i}" for i in range(n_cells)]


metadata = {'cell_type': np.random.choice(['A', 'B', 'C'], size=n_cells),
            'batch': np.random.choice(['1','2'], size=n_cells)}

metadata_df = pd.DataFrame(metadata, index=expression_df.index)

spatial_data = {'X': np.random.randint(0, 100, size=n_cells),
                'Y': np.random.randint(0, 100, size=n_cells)}
spatial_df = pd.DataFrame(spatial_data, index=expression_df.index)

image = np.random.randint(0, 256, size=(n_cells, 50, 50, 3)).astype(np.uint8)

# Create AnnData object
adata = create_anndata_from_dfs(expression_df, metadata_df, spatial_df, image)

if adata:
    print(adata)
    # Save the AnnData object
    adata.write_h5ad("merged_data.h5ad")
    print("AnnData object saved to merged_data.h5ad")

# Example without spatial data
adata_no_spatial = create_anndata_from_dfs(expression_df, metadata_df, image = image)
if adata_no_spatial:
    print(adata_no_spatial)
    adata_no_spatial.write_h5ad('no_spatial_data.h5ad')
    print("AnnData object with no spatial data saved to no_spatial_data.h5ad")

# Example without image data
adata_no_image = create_anndata_from_dfs(expression_df, metadata_df, spatial_df)
if adata_no_image:
    print(adata_no_image)
    adata_no_image.write_h5ad('no_image_data.h5ad')
    print("AnnData object with no image data saved to no_image_data.h5ad")

# Example with different image key
adata_diff_image_key = create_anndata_from_dfs(expression_df, metadata_df, spatial_df, image, image_key="my_image")
if adata_diff_image_key:
    print(adata_diff_image_key)
    adata_diff_image_key.write_h5ad('diff_image_key.h5ad')
    print("AnnData object with different image key saved to diff_image_key.h5ad")


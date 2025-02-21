import scanpy as sc
import os
import anndata
import matplotlib.pyplot as plt
import decoupler


def run(params, input_folder_paths, output_folder_paths, adata: anndata.AnnData):

    # parse online database using decoupler
    try:
        msigdb = decoupler.get_resource("MSigDB")
    except Exception as e:
        raise RuntimeError(f"Failed to fetch MSigDB resource: {e}")
    
    # Query GOBP for downstream analyses
    ## TODO: Parse collection from params
    metabo_database_name = params['metabo_database']
    metabo_database = msigdb[msigdb['collection'] == metabo_database_name]

    pathway_name = params['pathway_name']
    metabo_pathway = metabo_database[metabo_database['geneset'] == pathway_name]
    metabo_pathway = metabo_pathway[~metabo_pathway.duplicated()]

    # Select method
    method = params['method']
    # Run the selected method
    if method == "gsea":
        decoupler.run_gsea(
            adata,
            metabo_pathway,
            source="geneset",
            target="genesymbol",
            use_raw=False
        )
        # Extract the data from .obsm
        gsea_estimate = adata.obsm["gsea_estimate"]
        gsea_norm = adata.obsm["gsea_norm"]
        gsea_pvals = adata.obsm["gsea_pvals"]
        # Add UMAP coordinates as separate columns in .obs
        adata.obs["GSEA_score"] = gsea_estimate  # only support input one pathway
        adata.obs["GSEA_norm_score"] = gsea_norm
        adata.obs["GSEA_pvals"] = gsea_pvals


    elif method == "gsva":
        decoupler.run_gsva(
            adata,
            metabo_pathway,
            source="geneset",
            target="genesymbol",
            use_raw=False
        )
         # Extract the data from .obsm
        gsva_estimate = adata.obsm["gsva_estimate"]
        # Add UMAP coordinates as separate columns in .obs
        adata.obs["GEVA_score"] = gsva_estimate  # only support input one pathway


    elif method == "aucell":
        decoupler.run_aucell(
            adata,
            metabo_pathway,
            source="geneset",
            target="genesymbol",
            use_raw=False
        )
         # Extract the data from .obsm
        aucell_estimate = adata.obsm["aucell_estimate"]
        # Add UMAP coordinates as separate columns in .obs
        adata.obs["AUCell_score"] = aucell_estimate  # only support input one pathway


    elif method == "ora":
        decoupler.run_ora(
            adata,
            metabo_pathway,
            source="geneset",
            target="genesymbol",
            use_raw=False
        )
         # Extract the data from .obsm
        ora_estimate = adata.obsm["ora_estimate"]
        ora_pvals = adata.obsm["ora_pvals"]
        # Add UMAP coordinates as separate columns in .obs
        adata.obs["ORA_score"] = ora_estimate  # only support input one pathway
        adata.obs["ORA_pvals"] = ora_pvals  


    return adata
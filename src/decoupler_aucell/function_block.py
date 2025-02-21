import os
import scanpy as sc
import matplotlib.pyplot as plt
import decoupler


def run(params, input_folder_paths, output_folder_paths, adata):
    # parse online database using decoupler
    try:
        msigdb = decoupler.get_resource("MSigDB")
    except Exception as e:
        raise RuntimeError(f"Failed to fetch MSigDB resource: {e}")

    database = params['database']
    if database.lower() == "kegg":
        collection = 'kegg_pathways'
    elif database.lower() == "reactome":
        collection = 'reactome_pathways'
    elif database.lower() == "gobp":
        collection = 'go_biological_process'
    else:
        raise ValueError(f"Unsupported database: {database}")

    # Prepare computation parameters
    ## The element in pathway should be contained in geneset
    pathway = params['pathway']
    if isinstance(pathway, str):
        pathway = [pathway]

    # Query gene set from collection database for downstream analyses
    geneset = msigdb.query(f"collection == '{collection}'")
    geneset = geneset[geneset['geneset'].isin(pathway)]
    geneset = geneset[~geneset.duplicated()]

    decoupler.run_aucell(
        adata,
        geneset,
        source="geneset",
        target="genesymbol",
        use_raw=False,
    )

    # Export aucell result to adata format for visualization
    acts = decoupler.get_acts(adata, obsm_key='aucell_estimate')

    sc.pl.umap(acts, color=pathway, cmap='magma', vmax='p98', show=False)
    charts_path = output_folder_paths["charts"]
    umap_fig_path = os.path.join(charts_path, "umap.png")
    plt.tight_layout()
    plt.savefig(umap_fig_path, format='png')

    for p in pathway:
        adata.obs[p] = adata.obsm["aucell_estimate"][p]

    # Return the new AnnData object as computation result
    return adata
    
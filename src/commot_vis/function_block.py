import os
import matplotlib.pyplot as plt
import commot as ct


def run(params, input_folder_paths, output_folder_paths, adata):
    # parse params for querying database
    database = params['database']
    species = params['species']
    sig_type = params['signaling_type']

    plot_method = params['plot_method']
    background = params['background']
    cluster = params['clustering']

    ndsize = params['ndsize']
    scale = params['scale']
    normalize_v = params['normalize_v']
    normalize_v_quantile = params['normalize_v_quantile']
    grid_density = params['grid_density']

    # parse online database using decoupler
    try:
        lr_df = ct.pp.ligand_receptor_database(species=species, signaling_type=sig_type, database=database)
    except Exception as e:
        raise RuntimeError(f"Failed to fetch ligand receptorresource: {e}")

    # infer spatial communication
    ct.tl.spatial_communication(adata, database_name=database, df_ligrec=lr_df, dis_thr=500, heteromeric=True, pathway_sum=True)
    
    ct.tl.communication_direction(adata, database_name=database)
    ct.pl.plot_cell_communication(adata, database_name=database, plot_method=plot_method, background_legend=True,
        scale=scale, ndsize=ndsize, grid_density=grid_density, summary='sender', background=background, clustering=cluster, cmap='Alphabet',
        normalize_v=normalize_v, normalize_v_quantile=normalize_v_quantile)

    charts_path = output_folder_paths["charts"]
    umap_fig_path = os.path.join(charts_path, "ST_ccc.png")
    plt.tight_layout()
    plt.savefig(umap_fig_path, format='png')

    # Return the new AnnData object as computation result
    return adata
    
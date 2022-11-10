import yt
yt.enable_parallelism()

from yt.extensions.astro_analysis.halo_analysis import HaloCatalog

if __name__ == "__main__":
    es = yt.load_simulation('gas+dm-stars.enzo', 'Enzo')
    es.get_time_series()

    hc = HaloCatalog(
        data_ds=es, output_dir='rockstar_halos',
        finder_method='rockstar',
        finder_kwargs={'num_readers': 1, 'num_writers': 1})
    hc.create()


from setuptools import setup, find_packages


"""
pip install . -i https://mirrors.aliyun.com/pypi/simple --ignore-installed
"""

setup(
    name="magic_pipe",
    version="1.0.0",
    packages=find_packages(),
    author="Kai Li",
    python_requires='>=3.8',
    author_email="lemonsky79@gmail.com",
    install_requires=["colorama==0.4.4",
                      "tqdm==4.42.1",
                      "numpy==1.24.3",
                      "scipy==1.10.1",
                      "statsmodels==0.14.0",
                      "pygments==2.15.1",
                      "rpy2==3.4.5",
                      "sh==2.0.6",
                      "anndata==0.9.2",
                      'pandas==2.1.1'],
    entry_points={'console_scripts': ["magic=magic_pipe.magic:main"]},
    zip_safe=False,
    include_package_data=True,
    package_data={"magic_pipe": ["tools/scPBS_plot_cluster.R",
                                 "tools/scPBS_seurat_to_h5ad.R"]}
)

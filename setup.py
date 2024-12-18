from setuptools import setup, find_packages

setup(
    name="ivd_model",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "seaborn",
        "meshio",
        "SALib",
        "scikit-learn"
    ],
)
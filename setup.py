from setuptools import setup, find_packages

setup(
    name='msig',
    version='0.9.0',
    description='A small library to compute different per-vertex mesh signatures.',
    author="Dominik Penk",
    author_email="penk.dominik@gmail.com",
    packages=find_packages(include=['msig']),
    install_requires=[
        'numpy',
        'lapy',
        'trimesh',
        'scipy'
    ],
    url="https://github.com/DominikPenk/mesh-signatures",
    extras_require={
        'full': [
            'matplotlib',
            'imgui[glfw]',
            'imgui_datascience'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3.12.0',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    python_requires=">=3.12"
)
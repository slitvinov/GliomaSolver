from setuptools import Extension
from setuptools import setup

module = Extension(
    "mueler_brown",
    sources=["lib.cpp", "MRAG/MRAGcore/MRAGBoundaryBlockInfo.cpp"],
    extra_compile_args = ["-D_BLOCKSIZE_=8"],
)
setup(
    name="glioma_solver",
    version="1.0.0",
    description="Glioma solver module",
    url="https://github.com/slitvinov/GliomaSolver",
    author="Jana Lipková, Sergey Litvinov",
    author_email="slitvinov@gmail.com",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    ext_modules=[module],
)

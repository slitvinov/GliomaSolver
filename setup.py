from setuptools import Extension
from setuptools import setup

module = Extension(
    "mueler_brown",
    sources=["lib.cpp"],
)
setup(
    name="glioma_olver",
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

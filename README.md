## Input for the Reaction-Diffusion model
<dl>
  <dt>WM.dat, GM.dat</dt>
  <dd>white matter binary mask, gray matter binary mask</dd>
</dl>

## Model parameters
<dl>
  <dt>InputParameters.txt</dt>
  <dd> diffusivity in white matter [cm^2/day] (Dw) <br>
       proliferation rate [1/day] (rho) <br>
       final simulation time [day] (tend)
  </dd>

  <dt>TumorIC.txt</dt>
  <dd>tumor initial location (icx, icy, icz)</dd>
</dl>

## Output of the solver
<dl>
  <dt>HGG_data.dat</dt>
  <dd>tumor cell density [fraction]</dd>
</dl>

## Patient Data
<dl>
  <dt>tumFLAIR.dat</dt>
  <dd>binary mask, FLAIR MRI scan</dd>
</dl>

<dl>
  <dt>tumT1c.dat</dt>
  <dd>T1 gadolinium enhanced (T1Gd) scan, categorical 0, 1, 2, 4</dd>
</dl>

<dl>
  <dt>tumPET.dat</dt>
  <dd>PET-FET scan, float</dd>
</dl>


## Parameters of the likelihood
<dl>
  <dt>LikelihoodInput.txt</dt>
  <dd>PETsigma2, PETscale, T1uc, T2uc, slope</dd>
</dl>

## Install and run

```
$ python -m pip install -e .
```

## Run

```
$ ./fun.py
2.3475835715578407e+07
256x256x256le.raw
```

```
:; ./likelihood/likelihood
2.5588270414102585e+07
```

## Hacking

Coverage report
```
$ make -B 'CXXFLAGS = -fprofile-arcs -ftest-coverage' 'CFLAGS = -fprofile-arcs -ftest-coverage' 'LDFLAGS = -lgcov' 'LDDFLAGS = -lgcov'
$ ./brain
$ gcovr --html-details index.html
```

## Refernces

1. <http://tdo.sk/~janka/GliomaWebsite/index.html>

2. Lipková, J., Angelikopoulos, P., Wu, S., Alberts, E., Wiestler, B.,
   Diehl, C., ... & Menze, B. (2019). Personalized radiotherapy design
   for glioblastoma: Integrating mathematical tumor models, multimodal
   scans, and bayesian inference. IEEE transactions on medical
   imaging, 38(8), 1875-1884.
  [10.1109/TMI.2019.2902044](https://doi.org/10.1109/TMI.2019.2902044)

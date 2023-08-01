## Input for the Reaction-Diffusion model

<br>
<dl>
  <dt><CSF.dat, GM.dat, WM.dat</dt>
  <dd>white matter, gray matter,  cerebrospinal fluid</dd>
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
  <dd>tumor cell density?</dd>
</dl>

## Patient Data
<dl>
  <dt>tumFLAIR.dat, tumPET.dat, tumT1c.dat</dt>
  <dd>FLAIR MRI scan, PET-FET scan, T1Gd MRI scan</dd>
</dl>

## Parameters of the likelihood
<dl>
  <dt>LikelihoodInput.txt</dt>
  <dd>PETsigma2, PETscale, T1uc, T2uc, slope</dd>
</dl>

# Refernces

<http://tdo.sk/~janka/GliomaWebsite//index.html>

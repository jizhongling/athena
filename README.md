
Overview
--------
The reference detector at IP6 for Electron-Ion Collider experiment.

<a href="https://eicweb.phy.anl.gov/api/v4/projects/447/jobs/artifacts/master/raw/images/view01.pdf?job=report">
<img src="https://eicweb.phy.anl.gov/api/v4/projects/447/jobs/artifacts/master/raw/images/view01.png?job=report" width="400px" />
</a>

<br />
<a href="https://eicweb.phy.anl.gov/api/v4/projects/447/jobs/artifacts/master/raw/images/view01_top.pdf?job=report">
<img src="https://eicweb.phy.anl.gov/api/v4/projects/447/jobs/artifacts/master/raw/images/view01_top.png?job=report" width="400px" />
</a>


Getting Started
---------------

### Adding/changing detector geometry

Hint: **Use the CI/CD pipelines**.

To avoid dealing with setting up all the software, we recommend using the CI/CD to make changes.
Any feedback to help this process is appreciated.

Here is how to begin:

1. Look at existing detector constructions and reuse if possible. Note that "compact detector descriptions" -> xml files, and "detector construction" -> cpp file.
2. Modify xml file or detector construction. 
3. Create a WIP (or draft) merge request and look at the CI output for debugging. Then go to back to 2 if changes are needed.
4. Remove the WIP/Draft part of the merge request if you would like to see your changes merged into the master.

### Compiling 

First, see if the use case above is best for you. It most likely is and can save a lot of time for newcomers.
To run the simulation locally, we suggest using the singularity image.
More details can be found at the links below: 

- https://eic.phy.anl.gov/tutorials/eic_tutorial/
- https://eicweb.phy.anl.gov/containers/eic_container/






Related useful links
--------------------

- [Reference_Detector_doc] (https://escalate.readthedocs.io/projects/g4e/en/latest/detector_naming.html)
- [g4e] (https://gitlab.com/eic/escalate/g4e)
- [eic_tutorial] (https://argonne_eic.gitlab.io/tutorial/eic_tutorial/part2/adding_detectors/)
- [DD4hep] (https://github.com/AIDAsoft/DD4hep)
- [DD4hep_manual] (https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf)


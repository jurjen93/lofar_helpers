class: CommandLineTool
cwlVersion: v1.2
id: get_phasediff
label: Polarization Phase Difference
doc: This tool prepares measurement set for pulling phasediff scores from facetselfcal

baseCommand: python

inputs:
    - id: phasediff_ms
      type: Directory
      doc: Input measurement set
      inputBinding:
        position: 20
        valueFrom: $(inputs.phasediff_ms.path)
    - id: h5merger
      type: Directory
      doc: The h5merger directory.
    - id: selfcal
      type: Directory
      doc: The facetselfcal directory.

outputs:
    - id: phasediff_h5out
      type: File
      outputBinding:
        glob: "scalarphasediff*.h5"
    - id: logfile
      type: File[]
      outputBinding:
        glob: phasediff*.log

arguments:
  - valueFrom: $( inputs.selfcal.path + '/facetselfcal.py' )
  - valueFrom: -i phasediff
  - valueFrom: --forwidefield
  - valueFrom: --phaseupstations=core
  - valueFrom: --skipbackup
  - valueFrom: --uvmin=20000
  - valueFrom: --soltype-list=['scalarphasediff']
  - valueFrom: --solint-list=['10min']
  - valueFrom: --nchan-list=[6]
  - valueFrom: --docircular
  - valueFrom: --uvminscalarphasediff=0
  - valueFrom: --stop=1
  - valueFrom: --soltypecycles-list=[0]
  - valueFrom: --imsize=1600
  - valueFrom: --skymodelpointsource=1.0
  - valueFrom: --helperscriptspath=$(inputs.selfcal.path)
  - valueFrom: --helperscriptspathh5merge=$(inputs.h5merger.path)
  - valueFrom: --stopafterskysolve
  - valueFrom: --phasediff_only


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.phasediff_ms)
        writable: true
      - entry: $(inputs.h5merger)
      - entry: $(inputs.selfcal)

hints:
  - class: DockerRequirement
    dockerPull: vlbi-cwl
  - class: ResourceRequirement
    coresMin: 10


stdout: phasediff.log
stderr: phasediff_err.log

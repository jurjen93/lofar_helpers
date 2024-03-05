class: CommandLineTool
cwlVersion: v1.2
id: pre_averaging_dp3
label: DP3 Pre-averaging
doc: This step pre-averages measurement set for pulling phasediff scores from facetselfcal


baseCommand: DP3

inputs:
    - id: msin
      type: Directory
      doc: Input measurement set
      inputBinding:
        prefix: msin=
        position: 0
        valueFrom: $(inputs.msin.path)
        shellQuote: false
        separate: false

outputs:
    - id: phasediff_ms
      type: Directory
      outputBinding:
        glob: "*.phasediff.ms"
    - id: logfile
      type: File[]
      outputBinding:
        glob: dp3_prephasediff*.log

arguments:
  - steps=[avg]
  - msin.datacolumn=DATA
  - msout.storagemanager=dysco
  - avg.type=averager
  - avg.freqresolution=390.56kHz
  - avg.timeresolution=60
  - msout=$(inputs.msin.path+".phasediff.ms")

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.msin)
        writable: true

hints:
  - class: DockerRequirement
    dockerPull: vlbi-cwl
  - class: ResourceRequirement
    coresMin: 6


stdout: dp3_prephasediff.log
stderr: dp3_prephasediff_err.log

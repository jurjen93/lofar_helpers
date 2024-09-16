class: CommandLineTool
cwlVersion: v1.2
id: get_source_scores
label: Source Scores
doc: This tool determines the phasediff scores from h5 scalarphasediff input from facetselfcal

baseCommand: python

inputs:
    - id: phasediff_h5
      type: File[]
      doc: Phasedifference h5 from facetselfcal
      inputBinding:
        prefix: "--h5"
        position: 1
        itemSeparator: " "
        separate: true
    - id: selfcal
      type: Directory
      doc: The h5merger directory.


outputs:
    - id: phasediff_score_csv
      type: File
      outputBinding:
        glob: "*.csv"
    - id: logfile
      type: File[]
      outputBinding:
        glob: phasediff*.log


arguments:
  - valueFrom: $( inputs.selfcal.path + '/source_selection/phasediff_output.py' )


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.phasediff_h5)
        writable: true
      - entry: $(inputs.selfcal)

hints:
  - class: DockerRequirement
    dockerPull: vlbi-cwl


stdout: phasediff.log
stderr: phasediff_err.log

class: CommandLineTool
cwlVersion: v1.2
id: select_final
baseCommand: python3

inputs:
    - id: msin
      type: Directory[]
      doc: The input MS.
    - id: phasediff_score_csv
      type: File
      doc: The phase diff scores in a csv.

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.phasediff_score_csv)
      - entry: $(inputs.msin)
      - entryname: select_final.py
        entry: |
          import os
          import sys
          import pandas as pd
          from glob import glob
          import json

          inputs = json.loads(r"""$(inputs)""")
          msin = list(inputs['msin'])

          df = pd.read_csv(inputs['phasediff_score_csv')
          selection = df[df['spd_score'] < 2.4]['source'].to_list()

          os.system('mkdir -p output_ms')

          for s in selection:
            for ms in mslist:
              if s in ms:
                os.system(f"mv {ms} {ms}_selected.ms")

hints:
  - class: DockerRequirement
    dockerPull: vlbi-cwl

outputs:
    - id: selected_ms
      type: File
      outputBinding:
        glob: "*_selected.ms"
    - id: logfile
      type: File[]
      outputBinding:
        glob: select_final*.log

stdout: select_final.log
stderr: select_final.log

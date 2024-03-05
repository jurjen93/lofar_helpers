class: Workflow
cwlVersion: v1.2
id: selfcal pre-selection
doc: |
   This is a workflow to do a pre-selection for the LOFAR-VLBI pipeline direction-dependent calibrator selection

requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement

inputs:
    - id: msin
      type: Directory[]
      doc: The input MS.
    - id: h5merger
      type: Directory
      doc: The h5merger directory.
    - id: selfcal
      type: Directory
      doc: The selfcal directory.

steps:
    - id: dp3_prephasediff
      label: Pre-averaging with DP3
      in:
        - id: msin
          source: msin
      out:
        - id: phasediff_ms
      scatter: msin
      scatterMethod: dotproduct
      run: ./steps/dp3_prephasediff.cwl

    - id: get_phasediff
      label: Get phase difference with facetselfcal
      in:
        - id: phasediff_ms
          source: dp3_prephasediff/phasediff_ms
        - id: h5merger
          source: h5merger
        - id: selfcal
          source: selfcal
      out:
        - phasediff_h5out
      scatterMethod: dotproduct
      scatter: phasediff_ms
      run: ./steps/get_phasediff.cwl

    - id: get_selection_scores
      label: Calculate phase difference score
      in:
        - id: phasediff_h5
          source: get_phasediff/phasediff_h5out
        - id: h5merger
          source: h5merger
      out:
        - phasediff_score_csv
      run: ./steps/get_selection_scores.cwl

outputs:
    - id: phasediff_score_csv
      type: File
      outputSource: get_selection_scores/phasediff_score_csv

#!/bin/bash

WORKFLOW=$1

cwltool --pack $WORKFLOW > full_workflow.json
cwltool --print-dot full_workflow.json > workflow.dot
dot -Tsvg workflow.dot > $(basename ${WORKFLOW}).svg

rm full_workflow.json
rm workflow.dot

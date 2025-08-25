#!/bin/bash
source activate fme  
python -m fme.ace.inference inference_config${1}.yaml

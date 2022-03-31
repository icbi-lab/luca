#!/bin/bash

tar cvf - --dereference --hard-dereference containers | pigz -p 32 > data/zenodo/containers.tar.gz


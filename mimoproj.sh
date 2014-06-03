#!/bin/sh
exec /usr/sepp/bin/matlab -nojvm -nodisplay -nodesktop -nosplash -singleCompThread -r "$*"

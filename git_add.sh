#!/bin/bash

git add .
git reset -- Analysis/ContactStates/data4
git reset -- Analysis/ContactStates/data5
git reset -- Analysis/ContactStates/fig019
git reset -- Analysis/ContactStates/figInAc004less
git reset -- Analysis/ContactStates/figR05DefAgg
git reset -- Analysis/ContactStates/fig019
git reset -- Analysis/ContactStates/dataM5pN
git reset -- Analysis/ContactStates/figR05DefAgg_geq198
git reset -- Analysis/ContactStates/figR06DefAgg_geq054
git reset -- Analysis/ContactStates/figR07DefAgg_above019
git reset -- Analysis/ContactStates/figR07DefAgg_equal019
git reset -- Analysis/ContactStates/figR07DefAgg_geq019
git reset -- Analysis/ContactStates/figR08DefAgg_geq004
git reset -- Analysis/ContactStates/figR09DefAgg_geq001
git reset -- Analysis/ContactStates/figInAcAg
git reset -- Analysis/BindingMotors/data2
git status
echo "do:"
echo "git commit -m 'some comment' "
echo "git push -u origin master"
echo "Regards, Sirmaxford"

#===============================
# chmod -x git_add.sh
# source git_add.sh
# git commit -m "some comment"
# git push -u origin master
# Regards, sirmaxford
#===============================

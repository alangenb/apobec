#!/bin/bash

#A script to consolidate commands into a single task array for UGE submission.
#2015-09-04//Julian Hess
#
#Usage: qsubb.sh [commands files] [--pre=<prelude>] <UGE options>
#Commands can also (simultaneously) be given via STDIN.
#--pre can source a script to be run ahead of each task array job, e.g. for configuring environment.
#
#Examples: qsubb.sh commands1.txt commands2.txt -V -cwd -j y -l m_mem_free=4g -o log
#          echo -e "<commands>" | qsubb.sh <UGE options>
#          echo -e "<commands>" | qsubb.sh commands1.txt <UGE options>
#          qsubb.sh <commands> --pre=run_before.sh <UGE options>
#          qsubb.sh <commands> --pre=". /broad/software/scripts/useuse; reuse -q <dotkit>" <UGE options>

[ $# -eq 0 ] && { echo "Usage: qsubb.sh [command files] [--pre=<prelude script>] <UGE options>"; exit 1; }

#ensure we have UGER
type qsub | grep -q 'uge' || { echo "Error: UGE not on your path.  Exiting."; exit 1; }

#parse file arguments
while (($#)); do
	grep -q '^-' <<< $1 && break #stop at arguments
	if ! test -f "$1" && ! grep -q "^/dev/fd" <<< "$1"; then
		echo "Warning: not a regular file; skipping: $1";
		shift; continue
	fi
	cmdf[((f++))]=$1

	shift
done

#handle prelude
if grep -q '^--pre' <<< "$1"; then
	pre=`(sed -n 's/--pre="\?\(.*\)"\?/\1/p' | sed 's/[\]//gp') <<< $1`
	if test -f "$pre" || grep -q "^/dev/fd" <<< $pre; then
		pre=`cat "$pre"`
	fi
	shift
fi

#assemble commands
if [ ${#cmdf[@]} -gt 0 ]; then #read from files
	while read -r cat; do 
		cmd[((l++))]=$cat
	done < <(cat "${cmdf[@]}")
fi
if [ ! -t 0 ]; then
	while read -r cat; do
		cmd[((l++))]=$cat
	done < /dev/stdin
fi

[ ${#cmd[@]} -eq 0 ] && { echo "Nothing to do!"; exit 1; }

#verify that all commands exist before proceeding
for i in "${cmd[@]}"; do
	type `cut -d ' ' -f 1 <<< $i` > /dev/null 2>&1 || {
		echo "Error: command not found; dispatching halted: $i"
		exit 1
	}
done

#dispatch to qsub
echo 'cmd=`echo "'"`printf "%s\n" "${cmd[@]}"`"'" | sed -n ${SGE_TASK_ID}p`
echo '"$pre"'
echo ${SGE_TASK_ID}: $cmd
eval '"$pre"'
eval $cmd
echo Exit: $?' | qsub -t 1-${#cmd[@]} $@

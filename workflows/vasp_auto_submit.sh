#!/bin/bash
#===============================================================================
# Script Name: vasp_auto_submit.sh
# Purpose:
#   This script automates the submission of VASP jobs across multiple
#   calculation directories. It checks the status of each job by examining
#   the CONTCAR and OUTCAR files, and decides whether to resubmit or not.
#
#   For each directory:
#     1. If CONTCAR is missing or empty -> submit the job.
#     2. If CONTCAR exists and is non-empty:
#        a. If OUTCAR does NOT contain the success string -> backup POSCAR,
#           replace it with CONTCAR, and submit.
#        b. If OUTCAR contains the success string -> do nothing (job finished).
#
# Usage:
#   Modify the three configurable variables below (DIR_PATTERN, SUCCESS_STRING,
#   SUBMIT_CMD) to suit your needs, then execute the script:
#       ./vasp_auto_submit.sh
#
#   The script will loop over all directories matching DIR_PATTERN and process
#   each one as described.
#
#   NOTE: The script assumes that each directory contains POSCAR and OUTCAR
#         files, and that the submit command (e.g., ./sub.sh) is present in
#         each directory or available in PATH.
#===============================================================================

#------------------------------------------------------------------------------
# Configurable variables (modify as needed)
#------------------------------------------------------------------------------
# Pattern matching the VASP calculation directories (supports wildcards)
DIR_PATTERN="scan_x*"

# The success string to look for in OUTCAR (indicating a converged calculation)
SUCCESS_STRING="reached required accuracy - stopping structural energy minimisation"

# The command to submit a VASP job (must be executable from within the directory)
SUBMIT_CMD="./sub.sh"

#------------------------------------------------------------------------------
# Function: find_next_backup_number
#   Finds the next available two-digit number (01-99) for a backup file with
#   the given date stamp.
# Arguments:
#   $1 - date stamp in YYYYMMDD format
# Returns:
#   Echoes the next available number (two digits) and exits with status 0,
#   or exits with status 1 if all 01-99 are taken.
#------------------------------------------------------------------------------
find_next_backup_number() {
    local date_stamp="$1"
    local num
    for num in {01..99}; do
        local backup_file="POSCAR.bak.${date_stamp}${num}"
        if [ ! -f "$backup_file" ]; then
            echo "$num"
            return 0
        fi
    done
    echo "Error: All backup numbers (01-99) for date $date_stamp are used." >&2
    return 1
}

#------------------------------------------------------------------------------
# Main loop over directories
#------------------------------------------------------------------------------
# Enable globbing (wildcard expansion)
shopt -s nullglob

# Get current date for backup naming
CURRENT_DATE=$(date +%Y%m%d)

# Loop over all directories matching the pattern
for dir in $DIR_PATTERN; do
    # Ensure it is a directory
    [ -d "$dir" ] || continue

    echo "Processing directory: $dir"
    # Change into the directory (use subshell to avoid cd issues, but we need to perform actions)
    (
        cd "$dir" || { echo "Cannot enter $dir, skipping."; exit 1; }

        # Check if CONTCAR exists
        if [ ! -f "CONTCAR" ]; then
            echo "  CONTCAR not found. Submitting..."
            $SUBMIT_CMD
            exit 0
        fi

        # Check if CONTCAR is empty (size zero)
        if [ ! -s "CONTCAR" ]; then
            echo "  CONTCAR exists but is empty. Submitting..."
            $SUBMIT_CMD
            exit 0
        fi

        # CONTCAR exists and is non-empty: check OUTCAR
        # If OUTCAR does not exist, treat as not containing the success string
        if [ ! -f "OUTCAR" ] || ! grep -q "$SUCCESS_STRING" "OUTCAR"; then
            echo "  CONTCAR non-empty but OUTCAR missing or lacks success string."
            # Backup POSCAR with dated number
            # Find next available backup number
            next_num=$(find_next_backup_number "$CURRENT_DATE")
            if [ $? -ne 0 ]; then
                echo "  Failed to find backup number. Skipping this directory."
                exit 1
            fi
            backup_name="POSCAR.bak.${CURRENT_DATE}${next_num}"
            echo "  Renaming POSCAR to $backup_name"
            mv "POSCAR" "$backup_name"
            echo "  Replacing POSCAR with CONTCAR"
            mv "CONTCAR" "POSCAR"
            echo "  Submitting job..."
            $SUBMIT_CMD
        else
            echo "  CONTCAR non-empty and OUTCAR contains success string. Job is complete. No action."
        fi
    )
    # Note: the subshell ensures we return to the original directory automatically
done

echo "All directories processed."

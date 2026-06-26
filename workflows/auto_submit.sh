#!/bin/bash
#===============================================================================
# Script Name: auto_submit.sh
# Purpose:
#   This script automates the submission of DFT jobs across multiple
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
#   Additionally, after processing all directories, if *all* directories are
#   found to be finished (i.e., case 2b applies to every one), the script
#   generates a report file (default: completion_report.dat) listing each
#   directory name and the final E0 energy extracted from its OSZICAR file.
#   The report is only created when no submissions or restarts were performed.
#
# Usage:
#   Modify the configurable variables below (DIR_PATTERN, SUCCESS_STRING,
#   SUBMIT_CMD, REPORT_FILE) to suit your needs, then execute:
#       ./auto_submit.sh
#
#   The script will loop over all directories matching DIR_PATTERN and process
#   each one as described. After the loop, if all jobs are complete, it writes
#   the report in the current working directory.
#
#   NOTE: The script assumes that each directory contains POSCAR and OUTCAR
#         files, and that the submit command (e.g., ./sub.sh) is present in
#         each directory or available in PATH.
#===============================================================================

#------------------------------------------------------------------------------
# Configurable variables (modify as needed)
#------------------------------------------------------------------------------
# Pattern matching the calculation directories (supports wildcards)
DIR_PATTERN="scan_x*"

# The success string to look for in OUTCAR (indicating a converged calculation)
SUCCESS_STRING="reached required accuracy - stopping structural energy minimisation"

# The command to submit a job (must be executable from within the directory)
SUBMIT_CMD="./sub.sh"

# Name of the report file generated when all directories are finished.
# The file will contain two columns: directory name and the final E0 energy.
REPORT_FILE="completion_report.dat"

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

# Remove any stale incomplete marker from previous runs (if exists)
INCOMPLETE_MARKER=".auto_submit_incomplete"
rm -f "$INCOMPLETE_MARKER"

# Loop over all directories matching the pattern
for dir in $DIR_PATTERN; do
    # Ensure it is a directory
    [ -d "$dir" ] || continue

    echo "Processing directory: $dir"
    # Change into the directory (use subshell to avoid cd issues)
    (
        cd "$dir" || { echo "Cannot enter $dir, skipping."; exit 1; }

        # Check if CONTCAR exists
        if [ ! -f "CONTCAR" ]; then
            echo "  CONTCAR not found. Submitting..."
            $SUBMIT_CMD
            # Mark as incomplete (job submitted)
            touch "../$INCOMPLETE_MARKER"
            exit 0
        fi

        # Check if CONTCAR is empty (size zero)
        if [ ! -s "CONTCAR" ]; then
            echo "  CONTCAR exists but is empty. Submitting..."
            $SUBMIT_CMD
            touch "../$INCOMPLETE_MARKER"
            exit 0
        fi

        # CONTCAR exists and is non-empty: check OUTCAR
        # If OUTCAR does not exist, treat as not containing the success string
        if [ ! -f "OUTCAR" ] || ! grep -q "$SUCCESS_STRING" "OUTCAR"; then
            echo "  CONTCAR non-empty but OUTCAR missing or lacks success string."
            # Backup POSCAR with dated number
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
            touch "../$INCOMPLETE_MARKER"
        else
            echo "  CONTCAR non-empty and OUTCAR contains success string. Job is complete. No action."
            # No marker touched – this directory is finished
        fi
    )
    # Subshell ends – we are back in the original directory
done

#------------------------------------------------------------------------------
# Post-processing: generate report only if all directories are finished
#------------------------------------------------------------------------------
if [ ! -f "$INCOMPLETE_MARKER" ]; then
    echo "All directories appear to be finished. Generating energy report..."
    # Write header (optional) and two-column data
    {
        echo "# Directory   Final_E0 (from OSZICAR)"
        for dir in $DIR_PATTERN; do
            [ -d "$dir" ] || continue
            # Extract the last E0 value from OSZICAR if present
            if [ -f "$dir/OSZICAR" ]; then
                e0=$(grep "E0=" "$dir/OSZICAR" | tail -1 | awk '{print $NF}')
                if [ -z "$e0" ]; then
                    e0="NA"
                fi
            else
                e0="NA"
            fi
            echo "$dir   $e0"
        done
    } > "$REPORT_FILE"
    echo "Report written to: $REPORT_FILE"
else
    echo "Some directories are still incomplete (submissions or restarts were made)."
    echo "Report not generated. Remove '$INCOMPLETE_MARKER' manually if all jobs are actually done."
    # Clean up the marker after notifying
    rm -f "$INCOMPLETE_MARKER"
fi

echo "All directories processed."

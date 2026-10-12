#===============================================================================
# workflow_lib/common.sh -- functions shared by the numbered workflow drivers.
# Sourced, never executed. Drivers call wf_init first.
#
# Targets bash 4.2 (compy): set -Eeo pipefail (no -u: empty arrays), no
# namerefs, no `wait -n`. Never call a function that relies on `set -e` inside
# `if`, `||` or `&&`: bash then disables -e for everything it runs.
#===============================================================================

#===============================================================================
# logging
#===============================================================================
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%F %T')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

# Run a command, abort with a clear message if it fails.
must() { "$@" || die "command failed (rc=$?) in $(pwd): $*"; }

# ./xmlchange that aborts on failure (original scripts ignored xmlchange errors)
xmlc() { must ./xmlchange "$@" >/dev/null; }

# Source an environment script that may not be clean under `set -eu`.
source_env() {
  local f=$1 opts=$-
  [[ -f $f ]] || die "environment script not found: $f"
  set +eu
  # shellcheck disable=SC1090
  source "$f" >/dev/null
  [[ $opts == *e* ]] && set -e
  [[ $opts == *u* ]] && set -u
  return 0
}

#===============================================================================
# configuration
#===============================================================================
# wf_init <workflow dir>: source create_and_setup_case.sh, validate it, derive
# the settings used below and create runtmp/. Drivers pass the directory they
# resolved from SLURM_SUBMIT_DIR (batch) or their own location (interactive).
wf_init() {
  local dir=$1
  [[ -r ${dir}/create_and_setup_case.sh ]] || die "missing workflow configuration: ${dir}/create_and_setup_case.sh (submit from the workflow directory)"
  cd "$dir"
  # shellcheck disable=SC1091
  source ./create_and_setup_case.sh
  [[ $(readlink -f "$my_workflow_root") == "$(readlink -f "$dir")" ]] \
    || die "my_workflow_root (${my_workflow_root}) is not this directory (${dir})"
  check_submit_dir

  local v
  for v in my_workflow_root my_workflow_lib my_runtime_dir my_log_dir my_status_dir my_lock_dir \
           my_run_script_dir my_ensnum my_casedate my_casetod my_runtype my_compset my_resolution \
           my_machine my_project my_job_nnodes my_task_per_node my_jobqueue my_walltime my_e3sm_code \
           my_refcase my_refdir my_casename my_modeldir my_runpath my_stop_option my_stop_n; do
    [[ -n ${!v:-} ]] || die "create_and_setup_case.sh: ${v} is not set"
  done

  ENS_SIZE=${my_ensnum}       ENS_START=${my_enstart:-1}
  START_DATE=${my_casedate}   START_TOD=${my_casetod}
  RUN_TYPE=${my_runtype}      COMPSET=${my_compset}      RESOLUTION=${my_resolution}
  NODES=${my_job_nnodes}      TASKS_PER_NODE=${my_task_per_node}
  QUEUE=${my_jobqueue}        WALLTIME=${my_walltime}
  E3SM_CODE=${my_e3sm_code}
  REF_CASE=${my_refcase}      REF_DIR=${my_refdir%/archive}   # old configs pointed at .../archive
  CASE_NAME=${my_casename}    CASE_DIR=${my_modeldir}
  STOP_OPTION=${my_stop_option}  STOP_N=${my_stop_n}
  REST_OPTION=${my_rest_option:-nmonths}  REST_N=${my_rest_n:-1}
  MPAS_OCN_IC=${my_refocn_in:-} MPAS_ICE_IC=${my_refice_in:-}
  IC_SOURCE=${my_ic_source:-dart_analysis}   # dart_analysis | perturbed (Steps 0a/0b)
  LOG_DIR=${my_log_dir} STATUS_DIR=${my_status_dir} LOCK_DIR=${my_lock_dir} RUN_SCRIPT_DIR=${my_run_script_dir}
  RUN_TEMPLATE=${my_workflow_lib}/run_template/compile_and_setup_e3sm.${COMPSET}.sh

  [[ $ENS_SIZE =~ ^[0-9]+$ && $ENS_SIZE -ge 1 && $ENS_SIZE -le 99 ]] || die "my_ensnum must be 1..99"
  [[ $ENS_START =~ ^[0-9]+$ && $ENS_START -ge 1 && $ENS_START -le $ENS_SIZE ]] || die "my_enstart must be 1..my_ensnum"
  if ! [[ $START_DATE =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || ! date -d "$START_DATE" >/dev/null 2>&1; then
    die "my_casedate must be a valid YYYY-MM-DD date: ${START_DATE}"
  fi
  [[ $START_TOD =~ ^[0-9]{5}$ ]] && (( 10#$START_TOD < 86400 )) || die "my_casetod must be 5-digit seconds of day: ${START_TOD}"
  [[ $RUN_TYPE == AMIP || $RUN_TYPE == Full-CPL ]] || die "my_runtype must be AMIP or Full-CPL: ${RUN_TYPE}"
  [[ $IC_SOURCE == dart_analysis || $IC_SOURCE == perturbed ]] || die "my_ic_source must be dart_analysis or perturbed: ${IC_SOURCE}"
  [[ $STOP_OPTION =~ ^n(days|months|years)$ ]] || die "my_stop_option must be ndays, nmonths or nyears"
  [[ $STOP_N =~ ^[1-9][0-9]*$ ]] || die "my_stop_n must be a positive integer"
  [[ -r $RUN_TEMPLATE ]] || die "no run template for my_compset=${COMPSET}: ${RUN_TEMPLATE}"
  if [[ $RUN_TYPE == AMIP && $COMPSET == WCYCL* ]] || [[ $RUN_TYPE == Full-CPL && $COMPSET != WCYCL* ]]; then
    die "my_runtype=${RUN_TYPE} does not match my_compset=${COMPSET}"
  fi

  # time stamps used in restart file names / MPAS pointers
  START_HOUR=$(( 10#${START_TOD} / 3600 ))
  ATM_STAMP=${START_DATE}-${START_TOD}                       # eam/elm/mosart/cpl
  MPAS_STAMP=${START_DATE}_${START_TOD}                      # mpaso/mpassi
  MPAS_TIME=${START_DATE}_$(printf '%02d' "$START_HOUR"):00:00
  END_DATE=$(date -d "${START_DATE} + ${STOP_N} ${STOP_OPTION#n}" +%F)
  END_STAMP=${END_DATE}-${START_TOD}                         # end of the forecast

  mkdir -p "$LOG_DIR" "$STATUS_DIR" "$LOCK_DIR" "$RUN_SCRIPT_DIR"
}

# check_submit_dir: a batch job finds its workflow through SLURM_SUBMIT_DIR;
# refuse if the submitted script lives in a different workflow directory
# (e.g. `sbatch ../en10/2_run_e3sm_hindcast.sh` issued from another one).
check_submit_dir() {
  [[ -n ${SLURM_JOB_ID:-} ]] || return 0
  local cmd d
  cmd=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | head -1 || true)
  [[ $cmd == /* ]] || return 0
  d=$(cd "$(dirname "$cmd")" 2>/dev/null && pwd) || return 0
  d=$(readlink -f "$d")
  case $d in
    "$(readlink -f "$my_workflow_root")"|"$(readlink -f "$my_workflow_lib")"/*) ;;
    *) die "script ${cmd} belongs to another workflow than ${my_workflow_root}; cd to its directory and submit from there" ;;
  esac
}

summary_header() {
  log "workflow   ${my_workflow_root}"
  log "case       ${CASE_DIR}"
  log "reference  ${REF_DIR}  (${ATM_STAMP})"
  log "${RUN_TYPE} ${COMPSET}, ${ENS_SIZE} members, ${NODES} nodes x ${TASKS_PER_NODE} tasks each, queue ${QUEUE}, ${STOP_N} ${STOP_OPTION} (to ${END_STAMP})"
  if [[ $RUN_TYPE == Full-CPL ]]; then log "MPAS ICs   ${MPAS_OCN_IC%/*}"; fi
  log "members    $(tags "${MEMBERS[@]}")"
}

# normalize_bool <name> <value>: print TRUE/FALSE or die.
normalize_bool() {
  local v=${2^^}
  [[ $v == TRUE || $v == FALSE ]] || die "$1 must be TRUE or FALSE, got: $2"
  printf '%s' "$v"
}

tags() { local m; for m in "$@"; do printf '%s ' "$(mtag "$m")"; done; }
csv()  { local IFS=,; echo "$*"; }

#===============================================================================
# completion records (runtmp/status), key=value, written atomically
#===============================================================================
# Every record carries workflow_root, case, ensemble_size and valid_time; a
# record whose identity does not match this workflow is ignored, so records in
# a copied directory or from an earlier configuration can never be reused.
rec() { echo "${STATUS_DIR}/$1"; }

# record_write <name> [key=value ...]
record_write() {
  local f; f=$(rec "$1"); shift
  local tmp="${f}.tmp.$$"
  printf '%s\n' "valid_time=${ATM_STAMP}" "case=${CASE_NAME}" "ensemble_size=${ENS_SIZE}" \
    "workflow_root=${my_workflow_root}" "$@" "slurm_job_id=${SLURM_JOB_ID:-none}" \
    "written_at=$(date '+%F %T')" > "$tmp" || die "cannot write ${tmp}"
  mv -f "$tmp" "$f"
}

# record_get <name> <key>
record_get() { awk -F= -v k="$2" '$1 == k { sub(/^[^=]*=/, ""); print; exit }' "$(rec "$1")" 2>/dev/null || true; }

# record_valid <name>: true if the record exists and belongs to this workflow.
record_valid() {
  local f; f=$(rec "$1")
  [[ -s $f ]] || return 1
  [[ $(record_get "$1" workflow_root) == "$my_workflow_root" ]] || return 1
  [[ $(record_get "$1" case) == "$CASE_NAME" ]] || return 1
  [[ $(record_get "$1" ensemble_size) == "$ENS_SIZE" ]] || return 1
  [[ $(record_get "$1" valid_time) == "$ATM_STAMP" ]] || return 1
}

record_rm() { local n; for n in "$@"; do rm -f -- "$(rec "$n")"; done; }

# invalidate_after_run <m>: a new model run makes every downstream record stale.
invalidate_after_run() {
  local t; t=$(mtag "$1")
  record_rm "hindcast_complete.${t}" "hindcast_invalid.${t}" "compress_complete.${t}"
  rm -f -- "${STATUS_DIR}"/postprocess_complete.*."${t}".* 2>/dev/null || true
}

#===============================================================================
# case directory ownership: one workflow directory per case directory
#===============================================================================
# The case directory records which workflow directory created it, so a copied
# workflow that still points at the same case can never modify it.
owner_file() { echo "${CASE_DIR}/.hindcast_owner"; }

# case_owner: print the owning experiment directory ("" if none recorded).
case_owner() { cat "$(owner_file)" 2>/dev/null || true; }

# ownership_problem: print a problem description, or nothing if this
# experiment directory may manage CASE_DIR.
ownership_problem() {
  local owner; owner=$(case_owner)
  if [[ -n $owner && $owner != "$my_workflow_root" ]]; then
    echo "case ${CASE_DIR} belongs to experiment directory ${owner}; give this experiment its own case (my_refcase or my_casetag in create_and_setup_case.sh)"
  elif [[ -z $owner && -d $CASE_DIR ]]; then
    echo "case ${CASE_DIR} exists but was not created by this workflow (made by the old scripts?); to take it over set ADOPT_EXISTING_CASE=TRUE in 1_run_e3sm_ensemble_setup.sh"
  fi
}

claim_case_dir() {
  mkdir -p "$CASE_DIR"
  echo "$my_workflow_root" > "$(owner_file)"
}

# require_own_case: die if the case belongs to another experiment directory.
require_own_case() {
  local p; p=$(ownership_problem)
  if [[ $p == *"belongs to"* ]]; then die "$p"; fi
}

#===============================================================================
# members
#===============================================================================
mtag()       { printf 'EN%02d' "$1"; }
mdir()       { echo "${CASE_DIR}/$(mtag "$1")"; }
mcase()      { echo "${CASE_NAME}.$(mtag "$1")"; }
mscripts()   { echo "$(mdir "$1")/case_scripts"; }
mrundir()    { echo "$(mdir "$1")/run"; }
marchive()   { echo "$(mdir "$1")/archive"; }
mrefcase()   { echo "${REF_CASE}.$(mtag "$1")"; }

# select_members <spec>: fill the global array MEMBERS (member numbers).
#   all | 1,3,7-9 | EN03,EN05 | failed | unfinished
select_members() {
  local spec=${1:-all} part a b i st
  local -a picked=()
  case $spec in
    all) picked=($(seq "$ENS_START" "$ENS_SIZE")) ;;
    failed|unfinished)
      for i in $(seq "$ENS_START" "$ENS_SIZE"); do
        st=$(member_status "$i")
        case $spec:$st in
          failed:FAILED*|failed:INVALID|unfinished:FAILED*|unfinished:INVALID|unfinished:READY|unfinished:NO_CASE) picked+=("$i") ;;
        esac
      done ;;
    *)
      spec=${spec//EN/}
      for part in ${spec//,/ }; do
        if [[ $part =~ ^0*([0-9]+)-0*([0-9]+)$ ]]; then
          a=${BASH_REMATCH[1]} b=${BASH_REMATCH[2]}
        elif [[ $part =~ ^0*([0-9]+)$ ]]; then
          a=${BASH_REMATCH[1]} b=$a
        else
          die "bad member spec '${part}' (use all, failed, unfinished, or e.g. 1,3,7-9)"
        fi
        (( a >= 1 && b <= ENS_SIZE && a <= b )) || die "member range ${part} outside 1..${ENS_SIZE}"
        picked+=($(seq "$a" "$b"))
      done
      picked=($(printf '%s\n' "${picked[@]}" | sort -n -u)) ;;
  esac
  MEMBERS=("${picked[@]}")
}

#===============================================================================
# reference restarts
#===============================================================================
# ref_restdir <member>: directory holding the member's IC restarts, or empty.
ref_restdir() {
  local m=$1 d
  local probe; probe="$(mrefcase "$m").eam.i.${ATM_STAMP}.nc"
  for d in "${REF_DIR}/$(mtag "$m")/archive/rest/${ATM_STAMP}" \
           "${REF_DIR}/archive/rest/${ATM_STAMP}"; do
    if [[ -s ${d}/${probe} ]]; then echo "$d"; return 0; fi
  done
  return 1
}

# member_inputs <member> <restdir>: set atm_init lnd_init rof_init ice_init
# ocn_init and the list REQUIRED_INPUTS for this member.
member_inputs() {
  local m=$1 dir=$2 ref; ref=$(mrefcase "$m")
  atm_init=${dir}/${ref}.eam.i.${ATM_STAMP}.nc
  atm_rst=${dir}/${ref}.eam.r.${ATM_STAMP}.nc
  lnd_init=${dir}/${ref}.elm.r.${ATM_STAMP}.nc
  rof_init=${dir}/${ref}.mosart.r.${ATM_STAMP}.nc
  drv_rst=${dir}/${ref}.cpl.r.${ATM_STAMP}.nc
  # mpassi buildnml insists on <RUN_REFCASE>.mpassi.rst.* in RUNDIR for hybrid
  # runs, even when the streams point the model at a different file.
  ice_ref=${dir}/${ref}.mpassi.rst.${MPAS_STAMP}.nc
  if [[ $RUN_TYPE == Full-CPL ]]; then
    ice_init=${MPAS_ICE_IC}
    ocn_init=${MPAS_OCN_IC}
    REQUIRED_INPUTS=("$atm_init" "$lnd_init" "$rof_init" "$drv_rst" "$ice_ref" "$ice_init" "$ocn_init")
  else
    ice_init=$ice_ref
    ocn_init=
    REQUIRED_INPUTS=("$atm_init" "$lnd_init" "$rof_init" "$drv_rst" "$ice_init")
  fi
}

# check_member_inputs <member>: print problems (one per line), return 1 if any.
check_member_inputs() {
  local m=$1 dir f bad=0
  if ! dir=$(ref_restdir "$m"); then
    echo "$(mtag "$m"): no restart dir with $(mrefcase "$m").eam.i.${ATM_STAMP}.nc under ${REF_DIR}/{$(mtag "$m")/,}archive/rest/${ATM_STAMP}"
    return 1
  fi
  member_inputs "$m" "$dir"
  for f in "${REQUIRED_INPUTS[@]}"; do
    [[ -s $f && -r $f ]] || { echo "$(mtag "$m"): missing or empty: $f"; bad=1; }
  done
  # EAM starts from ncdata in a hybrid run; eam.r is only the rpointer.atm
  # target, so an empty placeholder (Step 0a) is enough.
  [[ -e $atm_rst ]] || { echo "$(mtag "$m"): missing: $atm_rst"; bad=1; }
  return $bad
}

#===============================================================================
# namelist / streams editing (idempotent)
#===============================================================================
# nl_set <user_nl file> <key> <value>: replace "key = ..." or append it.
nl_set() {
  local file=$1 key=$2 val=$3
  [[ -f $file ]] || die "namelist file not found: $(pwd)/${file}"
  if grep -qE "^[[:space:]]*${key}[[:space:]]*=" "$file"; then
    NL_KEY=$key NL_VAL=$val awk '
      BEGIN { k = ENVIRON["NL_KEY"]; v = ENVIRON["NL_VAL"] }
      $0 ~ "^[ \t]*" k "[ \t]*=" { print " " k " = " v; next }
      { print }' "$file" > "${file}.tmp" || die "failed to edit ${file}"
    mv "${file}.tmp" "$file"
  else
    printf ' %s = %s\n' "$key" "$val" >> "$file"
  fi
}

# stream_set <streams file> <stream name> <path>: point filename_template of the
# named stream at <path>. If the stream is missing and <path> is for the sea-ice
# restart_ic stream, the stream is added after the "restart" stream.
stream_set() {
  local file=$1 name=$2 path=$3
  [[ $path != *'"'* && $path != *'&'* ]] || die "unsupported characters in path: $path"
  if ! grep -q "name=\"${name}\"" "$file"; then
    [[ $name == restart_ic ]] || die "stream '${name}' not found in ${file}"
    ST_PATH=$path awk '
      { print }
      /<immutable_stream name="restart"/ { inb = 1 }
      inb && /\/>/ {
        inb = 0
        print ""
        print "<immutable_stream name=\"restart_ic\""
        print "                  type=\"input\""
        print "                  io_type=\"pnetcdf\""
        print "                  filename_template=\"" ENVIRON["ST_PATH"] "\""
        print "                  filename_interval=\"none\""
        print "                  input_interval=\"initial_only\" />"
      }' "$file" > "${file}.tmp" || die "failed to edit ${file}"
    mv "${file}.tmp" "$file"
  else
    ST_NAME=$name ST_PATH=$path awk '
      BEGIN { hit = 0 }
      $0 ~ "<(immutable_)?stream name=\"" ENVIRON["ST_NAME"] "\"" { inb = 1 }
      inb && /filename_template=/ {
        sub(/filename_template="[^"]*"/, "filename_template=\"" ENVIRON["ST_PATH"] "\""); hit = 1
      }
      { print }
      inb && /\/>/ { inb = 0 }
      END { if (!hit) exit 3 }' "$file" > "${file}.tmp" \
      || die "stream '${name}' in ${file} has no filename_template"
    mv "${file}.tmp" "$file"
  fi
  grep -qF "filename_template=\"${path}\"" "$file" || die "failed to set stream ${name} in ${file}"
}

#===============================================================================
# member status
#===============================================================================
# Cache of the user's active slurm jobs: "jobid jobname state" lines.
_SQUEUE_CACHE=
active_jobs() {
  if [[ -z $_SQUEUE_CACHE ]]; then
    _SQUEUE_CACHE=$(squeue -h -u "$USER" -o '%i %j %T' 2>/dev/null || true)
    [[ -n $_SQUEUE_CACHE ]] || _SQUEUE_CACHE=' '
  fi
  echo "$_SQUEUE_CACHE"
}

# member_status <member>: one of
#   NO_CASE   not set up              READY     set up, never submitted
#   PENDING   queued                  RUNNING   running
#   FINISHED  run ended OK, outputs not yet validated (3_check_hindcast.sh)
#   COMPLETE  outputs validated       INVALID   run ended but outputs incomplete
#   FAILED    run failed/cancelled    FAILED_ARCHIVE  run OK, st_archive failed
member_status() {
  local m=$1 t scripts name job last since
  t=$(mtag "$m"); scripts=$(mscripts "$m"); name=$(mcase "$m")
  if ! record_valid "setup_complete.${t}" || [[ ! -d $scripts ]]; then echo NO_CASE; return; fi
  job=$(active_jobs | awk -v n="$name" '$2 ~ "\\." n "$" || $2 == n { print $3; exit }' || true)
  case $job in
    PENDING|CONFIGURING|REQUEUED*|SUSPENDED) echo PENDING; return ;;
    ?*)                                      echo RUNNING; return ;;
  esac
  if record_valid "hindcast_complete.${t}"; then echo COMPLETE; return; fi
  if record_valid "hindcast_invalid.${t}"; then echo INVALID; return; fi
  if ! record_valid "hindcast_submitted.${t}"; then echo READY; return; fi
  # Only CaseStatus entries written after the last submission count.
  since=$(record_get "hindcast_submitted.${t}" submitted_at)
  last=$(awk -v t="$since" 'substr($0, 1, 19) >= t && / case\.(run|st_archive) (success|error|starting)/' \
           "${scripts}/CaseStatus" 2>/dev/null | tail -1 || true)
  case $last in
    *"case.st_archive success"*) echo FINISHED ;;
    *"case.st_archive error"*)   echo FAILED_ARCHIVE ;;
    *"case.run success"*)        echo FINISHED ;;
    *)                           echo FAILED ;;   # error, killed while running, or cancelled
  esac
}

# validate_member_outputs <m>: print problems with the archived forecast of a
# finished member (final restart + expected history files); return 1 if any.
validate_member_outputs() {
  local m=$1 archive name f d bad=0 stream comp
  archive=$(marchive "$m"); name=$(mcase "$m")
  for f in "${archive}/rest/${END_STAMP}/${name}.eam.r.${END_STAMP}.nc" \
           "${archive}/rest/${END_STAMP}/${name}.cpl.r.${END_STAMP}.nc"; do
    if [[ ! -s $f ]] || ! ncdump -h "$f" >/dev/null 2>&1; then echo "missing/unreadable final restart $f"; bad=1; fi
  done
  for stream in ${my_expected_daily_streams:-}; do
    comp=${stream%%:*}; stream=${stream#*:}
    for d in $(dates_in_range "$START_DATE" "$END_DATE"); do
      if ! ls "${archive}/${comp}/hist/${name}.${stream}.${d}-"*.nc >/dev/null 2>&1; then
        echo "missing ${stream} for ${d}"; bad=1; break
      fi
    done
  done
  for stream in ${my_expected_monthly_streams:-}; do
    comp=${stream%%:*}; stream=${stream#*:}
    for d in $(months_in_range "$START_DATE" "$END_DATE"); do
      f=${archive}/${comp}/hist/${name}.${stream}.${d}.nc
      if [[ ! -s $f ]] || ! ncdump -h "$f" >/dev/null 2>&1; then echo "missing/unreadable $f"; bad=1; fi
    done
  done
  return $bad
}

# dates_in_range <start> <end>: YYYY-MM-DD for each day in [start, end).
dates_in_range() {
  local d=$1
  while [[ $d < $2 ]]; do echo "$d"; d=$(date -d "$d + 1 day" +%F); done
}

# months_in_range <start> <end>: YYYY-MM for each month overlapping [start, end).
months_in_range() {
  local d; d=$(date -d "${1:0:7}-01" +%F)
  while [[ $d < $2 ]]; do echo "${d:0:7}"; d=$(date -d "$d + 1 month" +%F); done
}

#===============================================================================
# perturbed initial conditions (my_ic_source=perturbed, Steps 0a/0b)
#===============================================================================
# perturb_record_problem: print why the Step 0b record is unusable, or nothing.
perturb_record_problem() {
  local f=${STATUS_DIR}/perturb_complete.${ATM_STAMP}
  [[ -s $f ]] || { echo "no Step 0b record ${f}; run 0a_stage_initial_conditions.sh and 0b_perturb_initial_conditions.sh first"; return; }
  [[ $(awk -F= '$1 == "case" { print $2; exit }' "$f") == "$REF_CASE" ]] || { echo "Step 0b record ${f} is for another case"; return; }
  [[ $(awk -F= '$1 == "ensemble_size" { print $2; exit }' "$f") == "$ENS_SIZE" ]] || echo "Step 0b record ${f} has another ensemble size"
}

#===============================================================================
# preflight
#===============================================================================
# check_inputs: verify code, ownership and every selected member's initial
# conditions. Prints all problems, dies if there are any.
check_inputs() {
  local problems=() line m sst
  [[ -x ${E3SM_CODE}/cime/scripts/create_newcase ]] || problems+=("E3SM code not found: ${E3SM_CODE}")
  [[ -d $REF_DIR ]] || problems+=("reference case directory not found: ${REF_DIR}")
  if [[ $IC_SOURCE == perturbed ]]; then
    line=$(perturb_record_problem)
    if [[ -n $line ]]; then problems+=("$line"); fi
  fi
  for sst in $(sed -n 's/^readonly sst_\(data\|grid\)="\(.*\)".*/\2/p' "$RUN_TEMPLATE"); do
    [[ -f $sst ]] || problems+=("SST file named in $(basename "$RUN_TEMPLATE") not found: ${sst}")
  done
  if [[ -d $REF_DIR ]]; then
    for m in "${MEMBERS[@]}"; do
      while IFS= read -r line; do
        if [[ -n $line ]]; then problems+=("$line"); fi
      done < <(check_member_inputs "$m" || true)
    done
  fi
  line=$(ownership_problem)
  case $line in
    *"belongs to"*) problems+=("$line") ;;
    ?*)             warn "$line" ;;
  esac
  if (( ${#problems[@]} )); then
    for line in "${problems[@]}"; do printf '  PROBLEM: %s\n' "$line" >&2; done
    die "preflight failed with ${#problems[@]} problem(s)"
  fi
  log "preflight passed for $(tags "${MEMBERS[@]}")"
}

#===============================================================================
# cases: create, build, configure, verify, submit
#===============================================================================
# Layout (one directory per member):
#   ${CASE_DIR}/{case_scripts,build,run}           base case, builds e3sm.exe
#   ${CASE_DIR}/ENxx/{case_scripts,build,run,archive}

base_exe() { echo "${CASE_DIR}/build/e3sm.exe"; }

# materialize_run_template: copy the compset template to runtmp/run_scripts
# (the copy that is executed, kept for provenance).
materialize_run_template() {
  RUN_SCRIPT=${RUN_SCRIPT_DIR}/compile_and_setup_e3sm.sh
  cp -p -- "$RUN_TEMPLATE" "${RUN_SCRIPT}.tmp.$$" && chmod +x "${RUN_SCRIPT}.tmp.$$" \
    && mv -f "${RUN_SCRIPT}.tmp.$$" "$RUN_SCRIPT" || die "cannot materialize ${RUN_TEMPLATE}"
}

# run_case_script <case name> <root> <old exe or ""> : run the materialized
# compile_and_setup script for one case. Settings are passed in HC_* variables
# (no sed editing of the script). Needs member_inputs for the IC paths.
run_case_script() {
  local name=$1 root=$2 old_exe=$3
  [[ -n ${RUN_SCRIPT:-} && -x $RUN_SCRIPT ]] || die "run script not materialized"
  case $root in
    "$CASE_DIR"|"$CASE_DIR"/EN[0-9][0-9]) ;;
    *) die "refusing to set up a case outside ${CASE_DIR}: ${root}" ;;
  esac
  rm -rf -- "${root}/case_scripts"
  mkdir -p "$root"
  log "running $(basename "$RUN_SCRIPT") for ${name}"
  HC_MACHINE=$my_machine HC_PROJECT=$my_project HC_COMPSET=$COMPSET HC_RESOLUTION=$RESOLUTION \
  HC_CODE_ROOT=$E3SM_CODE HC_CASE_NAME=$name HC_CASE_ROOT=$root \
  HC_CASE_SCRIPTS_DIR=${root}/case_scripts HC_CASE_BUILD_DIR=${root}/build \
  HC_CASE_RUN_DIR=${root}/run HC_CASE_ARCHIVE_DIR=${root}/archive \
  HC_NODES=$NODES HC_TASK_PER_NODE=$TASKS_PER_NODE HC_JOB_QUEUE=$QUEUE HC_WALLTIME=$WALLTIME \
  HC_START_DATE=$START_DATE HC_START_TOD=$START_TOD \
  HC_STOP_OPTION=$STOP_OPTION HC_STOP_N=$STOP_N HC_REST_OPTION=$REST_OPTION HC_REST_N=$REST_N \
  HC_ATM_INIT=$atm_init HC_LND_INIT=$lnd_init HC_ROF_INIT=$rof_init HC_OCN_DATE=$MPAS_TIME \
  HC_OLD_EXECUTABLE=$old_exe \
    bash "$RUN_SCRIPT" || die "$(basename "$RUN_SCRIPT") failed for ${name}"
  [[ -x ${root}/case_scripts/xmlquery && -s ${root}/case_scripts/env_run.xml ]] \
    || die "incomplete case after setup: ${root}/case_scripts"
}

# build_base: create the base case and compile e3sm.exe (once per experiment).
build_base() {
  local dir
  dir=$(ref_restdir "$ENS_START") || die "no reference restarts for member ${ENS_START}"
  member_inputs "$ENS_START" "$dir"
  record_rm build_complete
  rm -rf -- "${CASE_DIR}/build"
  run_case_script "$CASE_NAME" "$CASE_DIR" ""
  [[ -x $(base_exe) ]] || die "build finished but $(base_exe) is missing"
  record_write build_complete "exe=$(base_exe)" "exe_mtime=$(stat -c %Y "$(base_exe)")"
  log "built $(base_exe)"
}

# setup_member <m>: create the member case with the base executable and
# configure it so it is ready to submit.
setup_member() {
  local m=$1 t root dir
  t=$(mtag "$m"); root=$(mdir "$m")
  record_valid build_complete && [[ -x $(base_exe) ]] || die "no completed build of $(base_exe)"
  dir=$(ref_restdir "$m") || die "${t}: no reference restarts"
  member_inputs "$m" "$dir"

  record_rm "setup_complete.${t}" "hindcast_submitted.${t}"
  invalidate_after_run "$m"
  run_case_script "$(mcase "$m")" "$root" "$(base_exe)"
  configure_member "$m"
  record_write "setup_complete.${t}" "member=${t}" "case_scripts=$(mscripts "$m")"
  log "${t}: case ready"
}

#-------------------------------------------------------------------------------
# configure_member <m>: point the member at its reference initial conditions.
# Idempotent; run before every submission.
#-------------------------------------------------------------------------------
configure_member() {
  local m=$1 dir f
  local scripts run ref
  scripts=$(mscripts "$m"); run=$(mrundir "$m"); ref=$(mrefcase "$m")

  dir=$(ref_restdir "$m") || die "$(mtag "$m"): no reference restarts found"
  member_inputs "$m" "$dir"
  for f in "${REQUIRED_INPUTS[@]}"; do
    [[ -s $f ]] || die "$(mtag "$m"): missing input $f"
  done
  [[ -e $atm_rst ]] || die "$(mtag "$m"): missing input $atm_rst"

  #--- case XML ---------------------------------------------------------------
  cd "$scripts"
  must ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=$(mcase "$m") \${EXEROOT}/e3sm.exe " >/dev/null
  xmlc RUN_TYPE=hybrid,CONTINUE_RUN=FALSE,RESUBMIT=0
  xmlc RUN_STARTDATE="$START_DATE",START_TOD="$START_TOD"
  xmlc STOP_OPTION="$STOP_OPTION",STOP_N="$STOP_N",REST_OPTION="$REST_OPTION",REST_N="$REST_N"
  xmlc GET_REFCASE=FALSE,RUN_REFCASE="$ref",RUN_REFDATE="$START_DATE",RUN_REFTOD="$START_TOD"
  xmlc RUN_REFDIR="$dir"
  xmlc DOUT_S=TRUE,DOUT_S_ROOT="$(marchive "$m")"
  xmlc JOB_QUEUE="$QUEUE"
  xmlc JOB_WALLCLOCK_TIME="$WALLTIME" --subgroup case.run

  #--- namelists ---------------------------------------------------------------
  local mpas_start; mpas_start="${START_DATE}_$(printf '%02d' "$START_HOUR")"
  nl_set user_nl_eam ncdata "\"${atm_init}\""
  nl_set user_nl_eam inithist "'MONTHLY'"
  nl_set user_nl_eam inithist_all .true.
  nl_set user_nl_elm finidat "\"${lnd_init}\""
  for f in check_finidat_year_consistency check_dynpft_consistency \
           check_finidat_fsurdat_consistency check_finidat_pct_consistency; do
    nl_set user_nl_elm "$f" .false.
  done
  nl_set user_nl_mosart finidat_rtm "\"${rof_init}\""
  nl_set user_nl_mpassi config_calendar_type "'gregorian'"
  nl_set user_nl_mpassi config_start_time "'${mpas_start}'"
  nl_set user_nl_mpassi config_initial_condition_type "'restart'"
  nl_set user_nl_mpassi config_do_restart .true.
  nl_set user_nl_mpassi config_restart_timestamp_name "'rpointer.ice'"
  if [[ $RUN_TYPE == Full-CPL ]]; then
    nl_set user_nl_mpaso config_calendar_type "'gregorian'"
    nl_set user_nl_mpaso config_do_restart .true.
    nl_set user_nl_mpaso config_start_time "'${mpas_start}'"
    nl_set user_nl_mpaso config_output_reference_time "'${MPAS_TIME}'"
    nl_set user_nl_mpaso config_restart_timestamp_name "'rpointer.ocn'"
  fi

  #--- run directory: restart pointers and links (before case.setup, see ice_ref)
  mkdir -p "$run"
  cd "$run"
  rm -f rpointer.*
  for f in "$atm_rst" "$lnd_init" "$rof_init" "$drv_rst" "$ice_ref"; do
    ln -sfn "$f" "$(basename "$f")"
  done
  if [[ $RUN_TYPE == Full-CPL ]]; then
    # mpaso buildnml wants <RUN_REFCASE>.mpaso.rst.*; the old scripts made it a
    # link to the MPAS ocean IC, which is also what streams.ocean points at
    ln -sfn "$ocn_init" "${ref}.mpaso.rst.${MPAS_STAMP}.nc"
  fi
  basename "$atm_rst"       > rpointer.atm
  echo "./$(basename "$lnd_init")" > rpointer.lnd
  echo "./$(basename "$rof_init")" > rpointer.rof
  basename "$drv_rst"       > rpointer.drv
  echo "$MPAS_TIME"         > rpointer.ice
  if [[ $RUN_TYPE == AMIP ]]; then
    # data ocean pointer, kept exactly as the original scripts wrote it
    local ocn1="${ref}.docn.r.${MPAS_STAMP}.nc" ocn2="${ref}.docn.rs1.${MPAS_STAMP}.bin"
    printf '%s\n%s\n' "$ocn1" "$ocn2" > rpointer.ocn
    for f in "$ocn1" "$ocn2"; do
      if [[ -f ${dir}/${f} ]]; then ln -sfn "${dir}/${f}" "$f"; fi
    done
  else
    echo "$MPAS_TIME" > rpointer.ocn
  fi

  cd "$scripts"
  must ./case.setup >/dev/null

  #--- MPAS streams (by stream name, kept in SourceMods) ------------------------------
  cd "$scripts"
  [[ -f ${run}/streams.seaice ]] || must ./preview_namelists >/dev/null
  patch_streams "$scripts" "$run" mpassi streams.seaice input "$ice_init" restart_ic "$ice_init"
  if [[ $RUN_TYPE == Full-CPL ]]; then
    patch_streams "$scripts" "$run" mpaso streams.ocean input "$ocn_init"
  fi

  must ./preview_namelists >/dev/null
  verify_member "$m"
}

# patch_streams <scripts> <run> <comp> <file> <stream> <path> [<stream> <path>...]
patch_streams() {
  local scripts=$1 run=$2 comp=$3 file=$4; shift 4
  local src=${scripts}/SourceMods/src.${comp}/${file}
  mkdir -p "$(dirname "$src")"
  if [[ ! -f $src ]]; then
    [[ -f ${run}/${file} ]] || die "no ${file} in ${run} to start from"
    cp "${run}/${file}" "$src"
  fi
  while (( $# >= 2 )); do
    stream_set "$src" "$1" "$2"
    shift 2
  done
  cp "$src" "${run}/${file}"
}

#-------------------------------------------------------------------------------
# verify_member <m>: confirm the generated namelists/streams/pointers really
# reference the intended initial conditions. Requires member_inputs.
#-------------------------------------------------------------------------------
verify_member() {
  local m=$1 run p tgt bad=0
  run=$(mrundir "$m")
  _expect() { grep -qF "$2" "${run}/$1" || { warn "$(mtag "$m"): ${1} does not reference $2"; bad=1; }; }
  _expect atm_in "$atm_init"
  _expect lnd_in "$lnd_init"
  _expect mosart_in "$rof_init"
  _expect streams.seaice "filename_template=\"${ice_init}\""
  if [[ $RUN_TYPE == Full-CPL ]]; then _expect streams.ocean "filename_template=\"${ocn_init}\""; fi
  for p in atm lnd rof drv; do
    tgt=$(head -1 "${run}/rpointer.${p}")
    [[ -e ${run}/${tgt} ]] || { warn "$(mtag "$m"): rpointer.${p} -> ${tgt} does not resolve"; bad=1; }
  done
  [[ -x $(mdir "$m")/build/e3sm.exe ]] || { warn "$(mtag "$m"): build/e3sm.exe missing"; bad=1; }
  (( bad == 0 )) || die "$(mtag "$m"): verification failed, not submitting"
}

# submit_member <m>: case.submit and record the submission.
submit_member() {
  local m=$1 t out jobid since
  t=$(mtag "$m")
  cd "$(mscripts "$m")"
  invalidate_after_run "$m"
  since=$(date '+%F %T')
  record_write "hindcast_submitted.${t}" "member=${t}" "submitted_at=${since}" "job_id=pending"
  set +e
  out=$(./case.submit 2>&1)
  local rc=$?
  set -e
  echo "$out"
  if (( rc != 0 )); then record_rm "hindcast_submitted.${t}"; die "${t}: case.submit failed (rc=${rc})"; fi
  jobid=$(echo "$out" | grep -oE 'Submitted job case\.run with id [0-9]+' | awk '{print $NF}' | tail -1 || true)
  record_write "hindcast_submitted.${t}" "member=${t}" "submitted_at=${since}" "job_id=${jobid:-unknown}"
  log "${t}: submitted case.run job ${jobid:-unknown}"
}

#===============================================================================
# parallel member processing
#===============================================================================
# run_pool <label> <function> <members...>: run <function> <member> for each
# member, at most JOBS at a time, each in its own subshell holding the member
# lock (runtmp/locks/member.ENxx.lock) and logging to runtmp/logs/<label>_ENxx.log.
# Sets POOL_OK / POOL_FAIL.
run_pool() {
  local label=$1 func=$2; shift 2
  local m rcdir logf
  rcdir=$(mktemp -d "${LOG_DIR}/.pool.XXXXXX")
  POOL_OK=() POOL_FAIL=()
  for m in "$@"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do sleep 2; done
    logf=${LOG_DIR}/${label}_$(mtag "$m").log
    log "$(mtag "$m"): ${label} started (log: ${logf})"
    {
      set +e   # keep going to record the status; the worker itself runs under -e
      (
        set -eo pipefail
        exec 9>"${LOCK_DIR}/member.$(mtag "$m").lock"
        flock -n 9 || die "$(mtag "$m") is locked by another workflow process"
        "$func" "$m"
      ) > "$logf" 2>&1
      echo $? > "${rcdir}/${m}"
    } &
  done
  wait
  for m in "$@"; do
    if [[ $(cat "${rcdir}/${m}" 2>/dev/null) == 0 ]]; then
      POOL_OK+=("$m")
    else
      POOL_FAIL+=("$m")
      warn "$(mtag "$m"): ${label} FAILED -- last lines of ${LOG_DIR}/${label}_$(mtag "$m").log:"
      tail -5 "${LOG_DIR}/${label}_$(mtag "$m").log" | sed 's/^/    /' >&2
    fi
  done
  rm -rf "$rcdir"
}

#===============================================================================
# driver helpers
#===============================================================================
# wf_exit_summary: EXIT trap for drivers; prints success/failure and run time.
WF_START_EPOCH=$(date +%s)
wf_exit_summary() {
  local rc=$? name; name=$(basename "${WF_SCRIPT_NAME:-$0}")
  if (( rc == 0 )); then
    echo "INFO: ${name} completed successfully in $(( $(date +%s) - WF_START_EPOCH ))s at $(date '+%F %T')"
  else
    echo "ERROR: ${name} failed with exit ${rc} after $(( $(date +%s) - WF_START_EPOCH ))s at $(date '+%F %T')" >&2
  fi
}

# stage_lock <name>: hold runtmp/locks/<name>.lock for the rest of the driver.
stage_lock() {
  exec 8>"${LOCK_DIR}/$1.lock"
  flock -n 8 || die "another $1 is already running (lock ${LOCK_DIR}/$1.lock)"
}

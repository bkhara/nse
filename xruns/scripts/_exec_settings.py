import numpy as np

# =======================================================================================
# case numbers
LidDrivenCavity = 201
FlowPastCylinder = 301

# =======================================================================================
# =======================================================================================
# program arguments
class ProgramArgumentStrings:
	PREP = "prep"
	RUN = "run"
	PROD = "prod"
	DEV = "dev"
	PP = "pp"
	CONS = "cons"
	CHECKDIV = "checkdiv"

# =======================================================================================
# =======================================================================================
# global file names
class ExecGlobalFileNames:
	Outfile = "out.txt"
	ExecLogFile = "multi_exec_log.txt"

# =======================================================================================
# =======================================================================================
# update this variable with the correct project path (i.e., /path/to/ns_vms)
PROJECT_PATH="${PROJECT_DIR}"

# =======================================================================================
# =======================================================================================
# multi options below
temporal_schemes = [
	"bdf",
	# "theta"
]

temporal_order = [
	# 1,
	2
]

class SupportingFiles:
	common = [
		"ns.petsc",
		"runex.sh",
		"job.sh",
		"job-dev.sh",
	]

	ldc = [
		"ldc_quad_bump.msh"
		# "generate_slice_configs_LDC.py"
	]

	fpc = [
		"fpc_uns.msh",
		"analyse_forces.py"
	]

class RunScript:
	simple = 'runex.sh'
	cluster_prod = 'job.sh'
	cluster_dev = 'job-dev.sh'

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

# =======================================================================================
# =======================================================================================
# MULTIRUN options below
class MMS_Conv_Data:
	rev = [
		# 1e0,
		# 1e1,
		# 1e2,
		1e3,
		1e4,
		1e5,
		1e6
	]
	nstepsv = [
		16,
		32,
		64,
		128,
		# 256,
	]
	bf = "quadratic"
	totalT = 2.5 * np.pi

class FPC2d_Data:
	rev = [
		1,
		15,
		30,
		60,
		100,
		150,
		200,
		250,
		300,
		400,
	]

class LDC2d_Data:
	rev = [
		100,
		400,
		1000,
		3200,
		5000,
		7500,
		10000
	]

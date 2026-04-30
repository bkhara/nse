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

	@staticmethod
	def decide_time_span_and_step(re):
		if re <= 50:
			totalT = 100
			dt = 0.2
		elif 50 < re <= 100:
			totalT = 150
			dt = 0.1
		elif 100 < re <= 200:
			totalT = 200
			dt = 0.05
		elif 200 < re <= 300:
			totalT = 250
			dt = 0.03
		elif 300 < re:
			totalT = 300
			dt = 0.01
		else:
			totalT = 300
			dt = 0.01
		return totalT, dt

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

	@staticmethod
	def decide_time_span_and_step(re):
		totalT = 50
		dt = 1
		if re < 5000:
			totalT = 5000 # 250 is sufficient
			dt = 1.0
		elif 5000 <= re < 7000:
			totalT = 5000 # 1500 is sufficient
			dt = 0.5
		elif 7000 <= re:
			totalT = 5000
			dt = 0.1
		return totalT, dt

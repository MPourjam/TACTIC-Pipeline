import shlex
from tactic_pipeline import processing_helper as proc_helper



PICRUST2_BIN = "/opt/conda/envs/picrust2/bin/picrust2_pipeline.py"

global PICRUSt2_logger
PICRUSt2_logger = proc_helper.gimmelogger(
    logger_name="run_tactic.picrust2"
)

PICRUSt2_logger.info("PICRUSt2 module loaded.")

def system_sub(*args, **kwargs):
    return proc_helper.system_sub(*args, logger_obj=PICRUSt2_logger, **kwargs)


def run_picrust2(*kwargs):
    """
    It receives the inact arguments passed to main CLI as string and adjust the input files if needed.
    Then it calls the main function of picrust2 module and returns its output.
    
    :param *kwargs: arguments passed to picrust2 main CLI as string
    """
    to_call = [
        PICRUSt2_BIN,
    ] + list(kwargs)

    PICRUSt2_logger.info(f"Calling PICRUSt2 with arguments: {' '.join(to_call)}")
    retcode = system_sub(to_call, force_log=True)
    if retcode != 0:
        PICRUSt2_logger.error(f"PICRUSt2 finished with non-zero exit code: {retcode}")
    else:
        PICRUSt2_logger.info("PICRUSt2 finished successfully.")
    return retcode

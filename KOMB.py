#!/usr/bin/env python
# See the LICENSE file included with this software for license information.

import os
import sys
import argparse
import logging
import time
from datetime import datetime
import subprocess
import shutil
import signal
import errno

__version__ = "2.0.0"

def handler(signum, frame):
    global SIGINT
    SIGINT = True
    sys.stderr.write("Caught request to terminate (CTRL+C), exiting...\n")
    sys.exit(128)

signal.signal(signal.SIGINT, handler)

class RunEnvironment():
    def __init__(self):
        self.CURRDATE = datetime.today().strftime('%Y-%m-%d')
        self.CURRTIME = datetime.today().strftime('%H:%M:%S')
        self.args     = self.parse_args()
        self.logger   = self.logging_setup()
        self.PATH_TO_KOMB = self.get_path_to_komb()
        
        if not self.PATH_TO_KOMB:
            sys.exit("Komb executable cannot be located, please add executable to $PATH or current directory.")

        try:
            os.mkdir(self.args.output_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and not self.args.overwrite:
                self.log(f"File or directory {self.args.output_dir} already exists.\nIf you want to overwrite use flag --overwrite.", 
                         logging.CRITICAL)
            elif exc.errno == errno.EEXIST:
                self.log(f"File or directory {self.args.output_dir} already exists and will be deleted.", 
                         logging.WARNING)
                if self.args.log_file not in ["stdout", "stderr"]:
                    os.remove(self.args.log_file)
                shutil.rmtree(self.args.output_dir)
                os.mkdir(self.args.output_dir)


    def get_path_to_komb(self):
        # First try to locate in PATH:
        try:
            result = subprocess.run(['which', 'komb2'], capture_output=True, text=True, check=True)
            komb2_path = result.stdout.strip()
            return komb2_path
        except subprocess.CalledProcessError:
            pass

        # Try to locate in current directory / src subdirectory:
        komb_directory = os.getcwd()
        komb2_path = os.path.join(komb_directory, 'komb2')

        if os.path.exists(komb_directory) and os.path.isfile(komb2_path):
            self.log(f"Using komb2 executable found at {komb2_path}. If this is the wrong one, " \
                            + "please add the correct executable to the $PATH", logging.WARNING)
            return komb2_path

        src_directory = os.path.join(komb_directory, 'src')
        komb2_path = os.path.join(src_directory, 'komb2')
        if os.path.exists(src_directory) and os.path.isfile(komb2_path):
            self.log(f"Using komb2 executable found at {komb2_path}. If this is the wrong one, " \
                            + "please add the correct executable to the $PATH", logging.WARNING)
            return komb2_path

        # Execute custom command to search for komb2 executable
        if os.name == 'posix':
            try:
                command = "find / -type f -name komb2 2>/dev/null"
                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                output, _ = process.communicate()
                output = output.decode().strip().split('\n')
                if output:
                    komb2_path = output[0]
                    self.log(f"Using komb2 executable found at {komb2_path}. If this is the wrong one, " \
                            + "please add the correct executable to the $PATH", logging.WARNING)
                    return komb2_path
                else:
                    self.log("Failed to find komb2 executable. Exiting...", logging.CRITICAL)
                    exit(-1)
            except FileNotFoundError:
                self.log("Failed to find komb2 executable. Exiting...", logging.CRITICAL)
                exit(-1)

        # Worst case scenario
        return None


    def setup(self):
        self.__init__()


    def parse_args(self):
        parser = argparse.ArgumentParser(description="""
        KOMB Analysis Pipeline example:
        python KOMB.py -i <read1.fq> -j <read2.fq> -k <k-mer size> -t <threads>
        """)

        io_args = parser.add_argument_group(title="Input/Output")
        io_args.add_argument(
            "-i",
            "--input-reads1",
            type=str,
            required=True,
            help="Path to the first sequencing reads file for paired-end data in FASTQ format"
        )
        io_args.add_argument(
            "-j",
            "--input-reads2",
            type=str,
            required=True,
            help="Path to the second sequencing reads file for paired-end data in FASTQ format"
        )
        io_args.add_argument(
            "-o",
            "--output-dir",
            type=str,
            default=f"komb_output_{self.CURRDATE}_{self.CURRTIME}",
            help="Path to the first sequencing reads file for paired-end data in FASTQ format"
        )
        io_args.add_argument(
            "--keep-alignments",
            action="store_true",
            help="Keep SAM files after the graph has been constructed (might require a lot of disk space)"
        )
        io_args.add_argument(
            "-e",
            "--log-file",
            type=str,
            default="stdout",
            help="File for logging [default: stdout]"
        )
        io_args.add_argument(
            "--overwrite",
            action="store_true",
            help="Delete <output-dir> and create a new one in case it exists"
        )

        common_args = parser.add_argument_group(title="Common arguments")
        common_args.add_argument(
            "-k",
            "--kmer-size",
            type=int,
            default=-1,
            required=True,
            help="k-mer size used for the *tig construction and subsequent analyses\nuse -1 to let KOMB automatically pick a value [default: -1]"
        )
        common_args.add_argument(
            "-t",
            "--num-threads",
            type=int,
            default=8,
            help="Maximum number of threads you want programs to use, note that some might use less than the amount specified"
        )
        common_args.add_argument(
            "-l",
            "--min-unitig-length",
            type=int,
            default=-1,
            help="Minimum length of a unitig to be kept for the analysis.\n \
                  Value -1 indicates setting this to the read length [default]\n \
                  Value 0 would result in keeping all unitigs, and values > 0 will apply the filter"
        )
        common_args.add_argument(
            "-v",
            "--verbosity",
            type=int,
            default=1,
            help="Logging level: 0 (DEBUG), 1 (INFO), 2 (ERROR)"
        )
        
        ggcat_args = parser.add_argument_group(title="GGCAT *tig construction")
        ggcat_args.add_argument(
            "-c",
            "--min-count",
            type=int,
            default=2,
            help="Minimum count required to keep a kmer [default: 2]"
        )
        ggcat_args.add_argument(
            "-m",
            "--ggcat-memory",
            type=int,
            default=8,
            help="Maximum memory usage for GGCAT (GB) [default: 8]"
        )
        ggcat_tig_type_args = parser.add_mutually_exclusive_group()
        ggcat_tig_type_args.add_argument(
            "--eulertigs",
            action="store_true",
            help="Generate Eulertigs instead of unitigs"
        )
        ggcat_tig_type_args.add_argument(
            "--greedy-matchtigs",
            action="store_true",
            help="Generate greedy matchtigs instead of unitigs"
        )
        ggcat_tig_type_args.add_argument(
            "--pathtigs",
            action="store_true",
            help="Generate pathtigs instead of unitigs"
        )

        bwa_args = parser.add_argument_group(title="BWA MEM parameters")
        bwa_args.add_argument(
            "--min-seed-length",
            type=int,
            default=18,
            help="Minimum seed length. Matches shorter than the value will be missed."
        )

        # komb_args = parser.add_argument_group(title="KOMB parameters")
        # komb_args.add_argument(

        # )
        return parser.parse_args()


    def logging_setup(self):
        logger = logging.getLogger("KOMB")
        logger.setLevel(logging.DEBUG)
        
        if self.args.log_file == "stdout":
            handler = logging.StreamHandler(stream=sys.stdout)
        elif self.args.log_file == "stderr":
            handler = logging.StreamHandler(stream=sys.stderr)
        else:
            handler = logging.FileHandler(self.args.log_file)
        
        error_handler = logging.StreamHandler(stream=sys.stderr)
        error_handler.setLevel(logging.ERROR)

        if self.args.verbosity < 1:
            handler.setLevel(logging.DEBUG)
        elif self.args.verbosity == 1:
            handler.setLevel(logging.INFO)
        else:
            handler.setLevel(logging.ERROR)

        formatter = logging.Formatter("%(asctime)s %(message)s",
                                      datefmt="%Y-%m-%d %I:%M:%S %p")
        handler.setFormatter(formatter)
        error_handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.addHandler(error_handler)
        return logger


    def log(self, msg, level=logging.DEBUG):
        if level == logging.DEBUG:
            self.logger.debug(msg)
        elif level == logging.INFO:
            self.logger.info(msg)
        elif level == logging.WARNING:
            self.logger.warning(msg)
        elif level == logging.ERROR:
            self.logger.error(msg)
        elif level == logging.CRITICAL:
            self.logger.error(msg)


    def RunGGCAT(self):
        cmd_str = f"ggcat build \
                    -k {self.args.kmer_size} \
                    -j {self.args.num_threads} \
                    -m {self.args.ggcat_memory} \
                    -s {self.args.min_count} \
                    -o {self.args.output_dir}/unitigs.fasta "
        if self.args.eulertigs:
            cmd_str += "--eulertigs "
        elif self.args.greedy_matchtigs:
            cmd_str += "--greedy-matchtigs "
        elif self.args.pathtigs:
            cmd_str += "--pathtigs "
        cmd_str += f"{self.args.input_reads1} {self.args.input_reads2}"

        p = subprocess.Popen(
            cmd_str,
            shell=True,
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            executable="/bin/bash",
            text=True
        )
        fstdout, fstderr = p.communicate()
        rc = p.returncode

        if rc != 0:
            self.logger.critical(f"GGCAT failed executing:\n{cmd_str}\n\nSTDOUT:{fstdout}\n\nSTDERR:{fstderr}")
            sys.exit(rc)
        else:
            self.log(fstdout, logging.DEBUG)
            self.log(fstderr, logging.DEBUG)


    def RunSeqKit(self):
        # NEED TO ADD COMPUTATION OF L FOR AUTO MODE
        cmd_str = f"seqkit seq \
                    -m {self.args.min_unitig_length} \
                    {self.args.output_dir}/unitigs.fasta"

        outfile = open(f"{self.args.output_dir}/unitigs.l{self.args.min_unitig_length}.fasta", "w")
        p = subprocess.Popen(
            cmd_str,
            shell=True,
            stdin=None,
            stdout=outfile,
            stderr=subprocess.PIPE,
            close_fds=True,
            executable="/bin/bash",
            text=True
        )
        fstdout, fstderr = p.communicate()
        rc = p.returncode

        if rc != 0:
            self.logger.critical(f"SeqKit failed executing:\n{cmd_str}\n\nSTDOUT:{fstdout}\n\nSTDERR:{fstderr}")
            sys.exit(rc)
        else:
            self.log(fstdout, logging.DEBUG)
            self.log(fstderr, logging.DEBUG)


    def RunBWAIndex(self):
        cmd_str = f"bwa index {self.args.output_dir}/unitigs.l{self.args.min_unitig_length}.fasta"

        p = subprocess.Popen(
            cmd_str,
            shell=True,
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            executable="/bin/bash",
            text=True
        )
        fstdout, fstderr = p.communicate()
        rc = p.returncode

        if rc != 0:
            self.logger.critical(f"BWA index failed executing:\n{cmd_str}\n\nSTDOUT:{fstdout}\n\nSTDERR:{fstderr}")
            sys.exit(rc)
        else:
            self.log(fstdout, logging.DEBUG)
            self.log(fstderr, logging.DEBUG)


    def RunBWAMem(self):
        def __bwamem(reads):
            cmd_str = f"bwa mem \
                        -a \
                        -k {self.args.min_seed_length} \
                        -t {self.args.num_threads} \
                        {self.args.output_dir}/unitigs.l{self.args.min_unitig_length}.fasta \
                        {reads}"

            outfile = open(f"{self.args.output_dir}/{reads.split('/')[-1]}.sam", "w")
            p = subprocess.Popen(
                cmd_str,
                shell=True,
                stdin=None,
                stdout=outfile,
                stderr=subprocess.PIPE,
                close_fds=True,
                text=True,
                executable="/bin/bash"
            )
            fstdout, fstderr = p.communicate()
            rc = p.returncode

            if rc != 0:
                self.logger.critical(f"BWA MEM failed executing:\n{cmd_str}\n\nSTDOUT:{fstdout}\n\nSTDERR:{fstderr}")
                sys.exit(rc)
            else:
                self.log(fstdout, logging.DEBUG)
                self.log(fstderr, logging.DEBUG)
        
        __bwamem(self.args.input_reads1)
        __bwamem(self.args.input_reads2)


    def RunKOMB(self):
        cmd_str = f"{self.PATH_TO_KOMB} \
                    -t {self.args.num_threads} \
                    -l {self.args.min_unitig_length} \
                    -o {self.args.output_dir} \
                    -i {self.args.output_dir}/{self.args.input_reads1.split('/')[-1]}.sam \
                    -j {self.args.output_dir}/{self.args.input_reads2.split('/')[-1]}.sam \
                    -u {self.args.output_dir}/unitigs.l{self.args.min_unitig_length}.fasta"

        p = subprocess.Popen(
            cmd_str,
            shell=True,
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            executable="/bin/bash",
            text=True
        )
        fstdout, fstderr = p.communicate()
        rc = p.returncode

        if rc != 0:
            self.logger.critical(f"KOMB-core failed executing:\n{cmd_str}\n\nSTDOUT:{fstdout}\n\nSTDERR:{fstderr}")
            sys.exit(rc)
        else:
            self.log(fstdout, logging.DEBUG)
            self.log(fstderr, logging.DEBUG)


    def RunCleanup(self):
        if not self.args.keep_alignments:
            os.remove(f"{self.args.output_dir}/{self.args.input_reads1.split('/')[-1]}.sam")
            os.remove(f"{self.args.output_dir}/{self.args.input_reads2.split('/')[-1]}.sam")


def main():
    run_env = RunEnvironment()
    
    t0 = time.time_ns()
    run_env.RunGGCAT()      # Run GGCAT to produce *tigs
    t1 = time.time_ns()
    run_env.log(f"Finished running GGCAT in {(t1-t0)/1000000.0:.3f} ms", logging.INFO)
    run_env.RunSeqKit()     # Only keep unitigs >= --min-unitig-length
    t2 = time.time_ns()
    run_env.RunBWAIndex()   # Build BWA index of the unitigs
    t3 = time.time_ns()
    run_env.RunBWAMem()     # Map reads back to the unitigs
    t4 = time.time_ns()
    run_env.log(f"Finished BWA index construction in {(t3-t2)/1000000.0:.3f} ms", logging.INFO)
    run_env.log(f"Finished read mapping in {(t4-t3)/1000000.0:.3f} ms", logging.INFO)
    t5 = time.time_ns()
    run_env.RunKOMB()       # Run KOMB constructing the graph and computing K-core decomposition
    t6 = time.time_ns()
    run_env.log(f"Finished running KOMB-core in {(t6-t5)/1000000.0:.3f} ms", logging.INFO)
    run_env.RunCleanup()    # If needed remove any files that are no longer used

    return 0


if __name__ == "__main__":
    main()
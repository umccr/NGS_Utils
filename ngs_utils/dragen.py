import csv
import glob
import json
import re
from collections import defaultdict
import yaml
from os import listdir
from os.path import join, abspath, pardir, splitext, basename, dirname, realpath, isdir, isfile, exists

from .file_utils import adjust_path, verify_dir, file_exists, safe_mkdir, verify_file, add_suffix
from .logger import critical, debug, info, err, warn
from .Sample import BaseSample, BaseBatch, BaseProject
from natsort import natsort_keygen


class DragenSample(BaseSample):
    natsort_key = natsort_keygen()

    def __init__(self, **kwargs):
        BaseSample.__init__(self, **kwargs)  # name, dirpath, work_dir, bam, vcf, phenotype, normal_match
        self.qc_files = []
        self.bam = None


class DragenBatch(BaseBatch):
    def __init__(
        self,
        tumor_normal_base_dir,
        tumor_normal_prefix,
        normal_base_dir,
        normal_prefix,
        **kwargs
    ):
        BaseBatch.__init__(self, **kwargs)
        self.somatic_caller = 'dragen'
        self.germline_caller = 'dragen'
        self.sv_caller = 'dragen'
        self.tumor_normal_base_dir = tumor_normal_base_dir
        self.normal_base_dir = normal_base_dir
        self.tumor_normal_prefix = tumor_normal_prefix
        self.normal_prefix = normal_prefix
        self.batch_qc_files = []

    def find_somatic_vcf(self, silent=False):
        somatic_vcf = str(self.tumor_normal_base_dir / f'{self.tumor_normal_prefix}.hard-filtered.vcf.gz')
        if isfile(somatic_vcf):
            verify_file(somatic_vcf, is_critical=True)
            if not silent:
                info(f'Found somatic VCF: {somatic_vcf}')
            self.somatic_vcf = somatic_vcf

    def find_germline_vcf(self, silent=False):
        germline_vcf = str(self.normal_base_dir / f'{self.normal_prefix}.hard-filtered.vcf.gz')
        if isfile(germline_vcf):
            verify_file(germline_vcf, is_critical=True)
            if not silent:
                info(f'Found germline VCF: {germline_vcf}')
            self.germline_vcf = germline_vcf

    def find_sv_vcf(self, silent=False):
        sv_vcf = str(self.tumor_normal_base_dir / f'{self.tumor_normal_prefix}.sv.vcf.gz')
        if isfile(sv_vcf):
            verify_file(sv_vcf, is_critical=True)
            if not silent:
                info(f'Found SV VCF: {sv_vcf}')
            self.sv_vcf = sv_vcf

    def find_qc_files(self):
        qc_files_shared = [
            '.vc_metrics.csv',
            '.ploidy_estimation_metrics.csv',
            '.fastqc_metrics.csv',
        ]
        qc_files_tumor = [
            *qc_files_shared,
            '.fragment_length_hist.csv',
            '.mapping_metrics.csv',
            '.wgs_contig_mean_cov_tumor.csv',
            '.wgs_coverage_metrics_tumor.csv',
            '.wgs_fine_hist_tumor.csv',
        ]
        qc_files_normal = [
            *qc_files_shared,
            '.wgs_fine_hist.csv',
            '.wgs_coverage_metrics.csv',
            '.wgs_contig_mean_cov.csv',
        ]
        # Construct full path of QC files and store
        for qc_file_suffix in qc_files_tumor:
            qc_file = self.tumor_normal_base_dir / f'{self.tumor_normal_prefix}{qc_file_suffix}'
            if not qc_file.exists():
                critical(f'DRAGEN QC file \'{qc_file}\' does not exist')
            self.tumors[0].qc_files.append(str(qc_file))
        for qc_file_suffix in qc_files_normal:
            qc_file = self.normal_base_dir / f'{self.normal_prefix}{qc_file_suffix}'
            if not qc_file.exists():
                critical(f'DRAGEN QC file \'{qc_file}\' does not exist')
            self.normals[0].qc_files.append(str(qc_file))

    def all_qc_files(self):
        return self.batch_qc_files + self.tumors[0].qc_files + self.normals[0].qc_files

    def add_tumor(self, prefix, rgid):
        sample = DragenSample(name=rgid, phenotype='tumor', batch=self, rgid=rgid)
        sample.bam = str(self.tumor_normal_base_dir / f'{prefix}_tumor.bam')
        self.tumors = [sample]
        if sample.name not in [s.name for s in self.parent_project.samples]:
            self.parent_project.samples.append(sample)
        return sample

    def add_normal(self, prefix, rgid):
        sample = DragenSample(name=rgid, phenotype='normal', batch=self, rgid=rgid)
        sample.bam = str(self.normal_base_dir / f'{prefix}.bam')
        self.normals = [sample]
        if sample.name not in [s.name for s in self.parent_project.samples]:
            self.parent_project.samples.append(sample)
        return sample


class DragenProject(BaseProject):

    def __init__(
        self,
        input_dirs,
        silent=False,
        include_samples=None,
        exclude_samples=None,
        genome_build=None,
        **kwargs
    ):
        BaseProject.__init__(self, **kwargs)
        self.genome_build = genome_build

        debug(f'Parsing project {", ".join(input_dirs)}')
        self.batch_by_name = DragenProject.find_batches(
            input_dirs,
            silent=silent,
            include_samples=include_samples,
            exclude_samples=exclude_samples,
            parent_project=self
        )

        if len(self.batch_by_name) == 1:
            self.project_name = list(self.batch_by_name.values())[0].name
        else:
            self.project_name = basename(input_dir)

    def add_batch(self, batch_name):
        batch = DragenBatch(name=batch_name, parent_project=self)
        self.batch_by_name[batch_name] = batch
        return batch

    @staticmethod
    def find_batches(runs, silent=False, include_samples=None, exclude_samples=None, parent_project=None):
        # Get identifiers and paths
        batch_name = runs['subject_id']
        rgsm_t = runs['tumor_normal_run']['tumor']
        tumor_normal_base_dir = runs['tumor_normal_run']['path']
        tumor_normal_prefix = runs['tumor_normal_run']['prefix']
        rgsm_n = runs['normal_run']['normal']
        normal_base_dir = runs['normal_run']['path']
        normal_prefix = runs['normal_run']['prefix']
        # Create batch with single sample
        batch = DragenBatch(
            tumor_normal_base_dir,
            tumor_normal_prefix,
            normal_base_dir,
            normal_prefix,
            name=batch_name,
            parent_project=parent_project,
        )
        # Add VCFs
        batch.find_somatic_vcf()
        batch.find_germline_vcf()
        batch.find_sv_vcf()
        # Add tumor/normal and normal samples with BAMs
        batch.add_tumor(tumor_normal_prefix, rgsm_t)
        batch.add_normal(normal_prefix, rgsm_n)
        # Add QC files
        batch.find_qc_files()
        # Return in format of 'batch_by_name'
        return {batch_name: batch}

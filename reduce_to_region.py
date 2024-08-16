import click
import os
import subprocess

@click.command()
@click.option('--region', help='Region to reduce to. Format chr:start-end', required=True)
@click.option('--gtf', help='Path to annotation GTF file', type=click.Path(exists=True), required=True)
@click.option('--output', help='Output directory', type=click.Path(exists=False), required=True)
@click.argument('bams', type=click.Path(exists=True), nargs=-1)
def main(region, gtf, bams, output):
    chr = region.split(':')[0]
    start, end = region.replace(',', '').split(':')[1].split('-')
    start, end = int(start), int(end)

    if not os.path.exists(output):
        os.makedirs(output)

    reduce_gtf(gtf, chr, start, end, output)
    for bam in bams:
        reduce_bam(bam, chr, start, end, output)

def reduce_gtf(gtf, chr, start, end, output):
    print(f'Reducing GTF file {gtf} to region {chr}:{start}-{end}')
    # reduce gtf file to region and output to output directory
    ext = os.path.splitext(gtf)[1]
    if ext == '.gz':
        cmd = f'zcat {gtf} | awk \'$1 == "{chr}" && $4 >= {start} && $5 <= {end}\' | gzip > {output}/reduced_{chr}_{start}_{end}.gtf.gz'
    elif ext == '.gtf':
        cmd = f'awk \'$1 == "{chr}" && $4 >= {start} && $5 <= {end}\' {gtf} > {output}/reduced_{chr}_{start}_{end}.gtf'
    print(cmd)
    subprocess.run(cmd, shell=True)

def reduce_bam(bam, chr, start, end, output):
    print(f'Reducing BAM file {bam} to region {chr}:{start}-{end}')
    assert os.path.exists(bam + '.bai'), 'BAM file must be indexed'
    # reduce bam file to region and output to output directory
    filename, ext = os.path.splitext(os.path.basename(bam))
    assert ext == '.bam', 'Input file must be a BAM file'
    cmd = f'samtools view -b {bam} {chr}:{start}-{end} > {output}/{filename}_{chr}_{start}_{end}.bam'
    print(cmd)
    subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    main()

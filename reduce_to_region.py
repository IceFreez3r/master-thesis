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
    output_bams = []
    for bam in bams:
        output_bams.append(reduce_bam(bam, chr, start, end, output))

    with open(os.path.join(output, f'metadata_{chr}_{start}_{end}.tsv'), 'w') as f:
        f.write('sample ID\tfile\n')
        for bam in output_bams:
            f.write(f'{'\t'.join(bam)}\n')

def reduce_gtf(gtf, chr, start, end, output):
    print(f'Reducing GTF file {gtf} to region {chr}:{start}-{end}')
    # reduce gtf file to region and output to output directory
    ext = os.path.splitext(gtf)[1]
    if ext == '.gz':
        cmd = f'zcat {gtf} | awk \'$1 == "{chr}" && $4 >= {start} && $5 <= {end}\' | gzip > {os.path.join(output, f'reduced_{chr}_{start}_{end}.gtf.gz')}'
    elif ext == '.gtf':
        cmd = f'awk \'$1 == "{chr}" && $4 >= {start} && $5 <= {end}\' {gtf} > {os.path.join(output, f'reduced_{chr}_{start}_{end}.gtf')}'
    print(cmd)
    subprocess.run(cmd, shell=True)

def reduce_bam(bam, chr, start, end, output):
    print(f'Reducing BAM file {bam} to region {chr}:{start}-{end}')
    assert os.path.exists(bam + '.bai'), 'BAM file must be indexed'
    # reduce bam file to region and output to output directory
    filename, ext = os.path.splitext(os.path.basename(bam))
    assert ext == '.bam', 'Input file must be a BAM file'
    output_file = os.path.join(output, f'{filename}_{chr}_{start}_{end}.bam')
    cmd = f'samtools view -b {bam} {chr}:{start}-{end} > {output_file}'
    print(cmd)
    subprocess.run(cmd, shell=True)
    cmd = f'samtools index {output_file}'
    print(cmd)
    subprocess.run(cmd, shell=True)
    return filename, os.path.abspath(output_file)

if __name__ == '__main__':
    main()

from astropy.io import fits
from glob import glob
import argparse


def update_fits_headers(imagelist, pipeline_repo, pipeline_commit, repro_repo):

    for image in imagelist:
        print('Updating FITS header:', image)

        with fits.open(image, mode='update') as hdul:
            header = hdul[0].header

            header.add_comment("======================================================================")
            header.add_comment(" The corresponding LOFAR data was processed and calibrated with PILOT:  ")
            header.add_comment(f"          {pipeline_repo} (commit {pipeline_commit})          ")
            header.add_comment("                                                                        ")
            header.add_comment("     For reproducibility for this particular dataset we refer to:       ")
            header.add_comment(f" {repro_repo}  ")
            header.add_comment("======================================================================")


def main():
    parser = argparse.ArgumentParser(description="Update FITS headers with provenance information")

    parser.add_argument("--fits", nargs="+", required=True, help="Input FITS files")
    parser.add_argument("--pipeline_repo", help="Pipeline GitHub URL",
                        default="https://github.com/LOFAR-VLBI/pilot")
    parser.add_argument("--pipeline_commit", help="Pipeline commit hash",
                        default="98ff43f")
    parser.add_argument("--repro_repo", help="Reproducibility repository URL",
                        default='https://github.com/LOFAR-VLBI/lofar_vlbi_helpers/tree/main/elais_200h')

    args = parser.parse_args()

    update_fits_headers(
        imagelist=args.fits,
        pipeline_repo=args.pipeline_repo,
        pipeline_commit=args.pipeline_commit,
        repro_repo=args.repro_repo
    )


if __name__ == "__main__":
    main()

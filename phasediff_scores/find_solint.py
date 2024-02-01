from source_selection.phasediff_output import GetSolint


if __name__ == "__main__":

    # set std score, for which you want to find the solint
    optimal_score = 2

    # reference solution interval in minutes
    ref_solint = 10

    # solution file
    h5 = '../P23872.h5'

    # get solution interval
    S = GetSolint(h5, optimal_score, ref_solint)
    solint = S.best_solint

    # OPTIONAL: plot fit
    S.plot_C("T=" + str(round(solint, 2)) + " min")

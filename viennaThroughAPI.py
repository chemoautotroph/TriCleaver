import RNA
from collections import defaultdict


def parse_dot_bracket(dotBracket):
    partners = dict()
    strands = dotBracket.split('&')
    num_strands = len(strands)
    strandNames = ['A', 'B']

    pairs = list()
    waitingToBePaired = list()
    waitingToBePairedPseudo = list()
    for strandId, strand in enumerate(strands):
        for baseId, baseBracketType in enumerate(strand):
            baseName = (strandId, baseId + 1)
            if baseBracketType == "(":
                waitingToBePaired.append(baseName)
            elif baseBracketType == "[":
                waitingToBePairedPseudo.append(baseName)
            else:
                if baseBracketType == ")":
                    partner = waitingToBePaired.pop()
                elif baseBracketType == "]":
                    partner = waitingToBePairedPseudo.pop()
                elif baseBracketType == ".":
                    partner = (None, None)

                pairs.append((partner,
                              baseName))  # non symmetric and out of order (for example, if we have (...), the () pair will be added after the dots
                partners[baseName] = partner
                if baseBracketType != ".":  # we don't want symmetry for the dummy (0,0) values
                    partners[partner] = baseName  # symmetry
    # print partners
    return partners


def fold(tasks, bool_300, folding_temp, num_strands):
    processed_results = list()
    for task in tasks:
        seq = task[0]
        constraint = task[1]

        # Create a model and set the temperature
        md = RNA.md()
        md.temperature = folding_temp

        # Prepare for cofolding if there are two strands
        if num_strands == 2:
            # Use cofolding
            fc = RNA.fold_compound(seq, md, RNA.OPTION_MFE | RNA.OPTION_PF | RNA.OPTION_HYBRID)
            fc.sc_add_bp(constraint)
            # (mfe_structure, mfe_energy) = fc.cofold()
            # RNA.update_cofold_params() Deprecated since version 2.6.3
            # (mfe_structure, mfe_energy) = fc.params_subst()
            fc.params_subst()  # Reset energy parameters to default values
            (mfe_structure, mfe_energy) = fc.mfe()
        else:
            # Use single strand folding
            fc = RNA.fold_compound(seq, md)
            fc.sc_add_bp(constraint)
            fc.params_subst()  # Reset energy parameters to default values
            (mfe_structure, mfe_energy) = fc.mfe()

        # Get base pair probabilities
        bppm = defaultdict(int)
        fc.exp_params_rescale(mfe_energy)
        fc.pf()
        for i in range(1, len(seq) + 1):
            for j in range(i + 1, len(seq) + 1):
                prob = fc.bpp(i, j)
                if prob > 0:
                    bppm[(i, j)] = prob
                    bppm[(j, i)] = prob

        # Calculate unpaired probabilities
        adj_dict = {i: [] for i in range(1, len(seq) + 1)}
        for (i, j), prob in bppm.items():
            adj_dict[i].append(j)
            adj_dict[j].append(i)

        for i in adj_dict.keys():
            paired_sum = sum(bppm[(i, j)] for j in adj_dict[i])
            unpaired_prob = 1 - paired_sum
            bppm[(i, None)] = unpaired_prob
            bppm[(None, i)] = unpaired_prob

        # Parse dot-bracket notation
        mfe_ad = parse_dot_bracket(mfe_structure)

        # Collect the results
        if num_strands == 1:
            ens_div = fc.mean_bp_distance()
            processed_results.append((mfe_structure, bppm, mfe_ad, mfe_energy, ens_div))
        else:
            processed_results.append((mfe_structure, bppm, mfe_ad, mfe_energy))

    return processed_results

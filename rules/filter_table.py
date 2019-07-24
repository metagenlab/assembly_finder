def select_assemblies(table, filter_by=False, nb=5, rank_to_select=None):
    if filter_by and rank_to_select != None:
        print('Filtering according to %s,%s and %s' % (rank_to_select, 'assembly status', 'category'))
        select_index = []
        unique_list = list(set(table[rank_to_select]))
        print(unique_list)

        if len(unique_list) > 1:
            for i in unique_list:
                select_index.append(table[table[rank_to_select] == i].sample(1).index[0])
        if len(unique_list) == 1:
            select_index = list(table.index)
    else:
        print('No filter specified, selecting according to assembly status and category')
        select_index = list(table.index)

    selection = table.loc[select_index]
    selection = selection.replace(
        ['reference genome', 'representative genome', 'Complete Genome', 'Chromosome', 'Contig',
         'Scaffold', 'na']
        , [0, 1, 2, 3, 4, 5, 6])
    sorted_sel = selection.sort_values(['Refseq_cat', 'AssemblySatus'], ascending=[True, True])[0:nb]
    return sorted_sel

'''
Main
'''

def get_plot(stackings):
    from matplotlib import pyplot as plt
    import mpld3
    fig, ax = plt.subplots()
    ax.scatter(
        [stacking["chi"] for stacking in stackings],
        [stacking["closest_distance"] for stacking in stackings],
    )
    ax.set_xlabel('Chi')
    ax.set_ylabel('Closest distance')
    return mpld3.fig_to_html(fig)

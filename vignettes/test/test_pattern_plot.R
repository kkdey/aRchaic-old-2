

pattern_plot(file="../summary_data/I0026.q30.csv",
             pattern="C->T",
             plot_type="left")

pattern_plot(file="../summary_data/I0026.q30.csv",
             pattern="G->A",
             plot_type="right")

pattern_plot(file="../summary_data/I0026.q30.csv",
             pattern="T->AAG",
             plot_type="both")

pattern_plot(file="../summary_data/I0026.q30.csv",
             pattern="T->A",
             plot_type="right")


# test example
##
out <- damage.build.counts(file="../summary_data/87_I_RISE_study/I0026.q30.csv")

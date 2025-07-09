# extended version of logo_shared_gkmer.py; used to generate Fig 5C.
import logomaker
import numpy as np
import pandas as pd
import sys
from matplotlib import pyplot as plt
def main(argv = sys.argv):
    
    if(len(argv)!= 9):
        print("Usage <0> <matrix file 1> <dsvm 1> <ID for matrix file 1> <matrix file 2> <dsvm 2> <ID for matrix file 2>  <info> <ofile name>")
        sys.exit(0)
    fig, (ax_q, ax_t) = plt.subplots(nrows=2)

    mfile1 = argv[1]
    dsvm1 = open(argv[2], 'r')
    id1 = argv[3]
    mfile2 = argv[4]
    dsvm2 = open(argv[5], 'r')
    id2 = argv[6]
    cons_info = argv[7]
    ofname = argv[-1]


    # conservation logo track
    mat_q = pd.read_csv(mfile1, delim_whitespace=True, index_col=0, comment='#')
    logo_q = logomaker.Logo(mat_q, ax = ax_q, center_values=False, figsize = [100, 5])

    mat_t = pd.read_csv(mfile2, delim_whitespace=True, index_col=0, comment='#')
    logo_t = logomaker.Logo(mat_t, ax = ax_t, center_values=False, figsize = [100, 5])

    ymin, ymax = ax_q.get_ylim()
    xmin, xmax = ax_q.get_xlim()
    xpos = (xmax-xmin)*0.01 + xmin
    ypos = (ymax-ymin)*1.4 + ymin
    ax_q.text(xpos, ypos, cons_info, fontsize=25)
    ax_q.set_title(id1, fontsize=25)
    ax_t.set_title(id2, fontsize=25)

    ax_q.set_xlabel("position", fontsize=20)
    ax_t.set_xlabel("position", fontsize=20)

    ax_q.set_ylabel("gkm cons.", fontsize=20)
    ax_t.set_ylabel("gkm cons.", fontsize=20)

    ax_q.set_ylim([-((mat_q.values).max()*2), (mat_q.values).max()*2])
    ax_t.set_ylim([-((mat_t.values).max()*2), (mat_t.values).max()*2])

    # also add gkmsvm score track 
    ax_q_svm = ax_q.twinx();
    models_q = [] # [(model_name, vector of score)]
    abs_max_q = -100 #useful for plotitng later
    for i, line in enumerate(dsvm1):
        words = line.split()
        if(i==0):
            for model_id in words[1:]:
                models_q.append([model_id, []])
        else:
            for j, sc in enumerate(words[1:]):
                models_q[j][1].append(eval(sc))
                if(abs(eval(sc)) > abs_max_q):
                    abs_max_q = abs(eval(sc))
    for i in range(0,len(models_q)):
        x = np.array(range(0,len(models_q[i][1])))
        y = models_q[i][1]  
        ax_q_svm.plot(x, y, linewidth = 2, label=models_q[i][0])# label = models_q[i][0])
        ax_q_svm.fill_between(x, y, alpha = 0.2)
    ax_q_svm.legend(loc='upper right', fontsize=15)
    ax_q_svm.set_ylabel("avg delta-SVM", fontsize=20)
    ax_q_svm.set_ylim([-abs_max_q*2, abs_max_q*2])

    # same for target 
    ax_t_svm = ax_t.twinx();
    models_t = [] # [(model_name, vector of score)]
    abs_max_t = -100
    for i, line in enumerate(dsvm2):
        words = line.split()
        if(i==0):
            for model_id in words[1:]:
                models_t.append([model_id, []])
        else:
            for j, sc in enumerate(words[1:]):
                models_t[j][1].append(eval(sc))
                if(abs(eval(sc)) > abs_max_t):
                    abs_max_t = abs(eval(sc))
    for i in range(0, len(models_t)):
        x = np.array(range(0,len(models_t[i][1])))
        y = models_t[i][1] 
        ax_t_svm.plot(x, y, linewidth = 2, label = models_t[i][0])
        ax_t_svm.fill_between(x, y, alpha = 0.2)
    ax_t_svm.legend(loc='upper right', fontsize=15)
    ax_t_svm.set_ylabel("avg delta-SVM", fontsize=20)
    ax_t_svm.set_ylim([-abs_max_t*2, abs_max_t*2])



    fig.set_size_inches(30, 12)
    fig.tight_layout()
    fig.savefig(ofname)

main()

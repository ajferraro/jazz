import matplotlib.pyplot as plt

def plot_styling(mode='jazz_paper'):
    try:
        plt.style.use(mode)
    except:
        print 'No style sheet available'


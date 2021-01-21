from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
import argparse


def venn(vennset, vennlabel, title, venncolor=('#990000','#00909e', '#ff971d'), vennalpha=0.5):
	fig = plt.figure(figsize=(4,4))
	if len(vennset) == 7:
		venn3(subsets=vennset, set_labels=vennlabel, set_colors=venncolor, alpha=vennalpha)
		# Save figure 
		plt.savefig('venn3.png', format='png', bbox_inches='tight', dpi=300)
		# Generate the Markdown preview
        #print('![picture](venn3.png)')
	elif len(vennset) == 3:
		v = venn2(subsets=vennset, set_labels=vennlabel, set_colors=venncolor, alpha=vennalpha)
		# Save figure 
		v.get_label_by_id('11').set_text('%0.0000f%%' % (vennset[2]/sum(vennset)*100))
		v.get_label_by_id('01').set_text('%0.0000f%%' % (vennset[1]/sum(vennset)*100))
		v.get_label_by_id('10').set_text('%0.0000f%%' % (vennset[0]/sum(vennset)*100))
		plt.title(title)
		plt.savefig('venn2.png', format='png', bbox_inches='tight', dpi=400)


        # Generate the Markdown preview
        #print('![picture](venn2.png)')
	else:
		print("Error: check the set dataset")

# Load input data
parser = argparse.ArgumentParser()
parser.add_argument('--set', help = 'Set of data for Venn diagram: order the set in "100,010,110,001,101,011,111" format', dest = 'set')
parser.add_argument('--labels', help = ' Labels for each of the groups', dest = 'labels')
parser.add_argument('--title', help = 'Title for the plot', dest = 'title')
args = parser.parse_args()

vennset = args.set.split(',')
vennset  = list(map(int, vennset))
labels = args.labels.split(',')
title = args.title
venn(vennset,labels, title )


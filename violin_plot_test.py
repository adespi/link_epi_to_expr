for x in range(-5,6):
   if len(Y[np.round(X)==x])>0:
      a = plt.violinplot(Y[np.round(X)==x],[x])
      for pc in a['bodies']:
         pc.set_color('C0')
      a['cbars'].set_color('C0')
      a['cmins'].set_color('C0')
      a['cmaxes'].set_color('C0')

plt.show()

plt.violinplot([Y[np.round(X)==x]for x in range(-2,3)],range(-2,3)); plt.show()

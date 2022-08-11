import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.add_feature(cfeature.LAND, alpha = 0)
ax.add_feature(cfeature.OCEAN, color = 'White')
ax.coastlines(color = 'Black',resolution='110m')
ax.set_global()
plt.tight_layout()
plt.show()
fig.savefig('/Users/chrishancock/Desktop/temp.png', transparent=True,dpi=600,bbox_inches='tight')

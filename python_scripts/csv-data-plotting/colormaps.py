import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def save_middle_range_cmap(cmap_name, start=0, end=1, steps=10, filename=None):
    # Get the colormap
    cmap = cm.get_cmap(cmap_name)

    # Get a range of values around the middle of the colormap
    colors = cmap(np.linspace(start, end, steps))

    # Reshape to make it a 2D array for visualization
    colors_image = np.vstack([colors] * steps)  # Repeat the colors to create an image

    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 1))

    # Display the image
    ax.imshow(colors_image, aspect="auto")

    # Hide the axes
    ax.axis('off')

    # Default filename if not provided
    if filename is None:
        filename = f'{cmap_name}_middle_range.png'

    # Save as PNG
    plt.savefig('colormaps/' + filename, format='png', bbox_inches='tight')
    plt.close()

    print(f"Image has been saved to colormaps/{filename}")


def visualize_cmap(cmap_name, filename=None):
    # Get the colormap
    cmap = cm.get_cmap(cmap_name)

    # Create an array representing the full range of the colormap
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    gradient = np.vstack((gradient, gradient))

    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 1))

    # Display the gradient as an image
    ax.imshow(gradient, aspect='auto', cmap=cmap)

    # Hide the axes
    ax.set_axis_off()

    # Default filename if not provided
    if filename is None:
        filename = f'{cmap_name}_full_range.png'

    # Save as PNG
    plt.savefig('colormaps/' + filename, format='png', bbox_inches='tight')
    plt.close()

    print(f"Image has been saved to colormaps/{filename}")


colormaps = [
    'viridis', 'inferno', 'plasma', 'magma',
    'RdBu', 'PiYG', 'coolwarm',
    'twilight', 'hsv',
    'tab10', 'tab20', 'Set1',
    'rainbow', 'jet'
]

fig, axes = plt.subplots(1, len(colormaps), figsize=(20, 2))
for ax, cmap in zip(axes, colormaps):
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    gradient = np.vstack((gradient, gradient))
    ax.imshow(gradient, aspect="auto", cmap=plt.get_cmap(cmap))
    ax.set_title(cmap)
    ax.axis('off')

# Save the figure as a PNG file
png_filename = 'colormaps.png'
plt.savefig(png_filename, format='png', bbox_inches='tight')
plt.close()

# query middle range of colormaps
save_middle_range_cmap('inferno')
save_middle_range_cmap('RdBu')


# Example usage
visualize_cmap('inferno')


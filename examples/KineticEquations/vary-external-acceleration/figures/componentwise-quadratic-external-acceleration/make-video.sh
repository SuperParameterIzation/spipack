ffmpeg -r 50 -i frames/output-%12d.png -vcodec libx264 -y -an simulation.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"

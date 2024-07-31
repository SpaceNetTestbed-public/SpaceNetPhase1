timestamp=$(date +%s)
echo $timestamp
wget -O "./utils/starlink_tles/starlink_${timestamp}" https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle


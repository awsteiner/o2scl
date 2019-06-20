qput:
	rsync -Cavzun Hawaii2017/ awsteiner@lit:Hawaii2017
	rsync -Cavzun Hawaii2017Edit/ awsteiner@lit:Hawaii2017Edit

put:
	rsync -Cavzu Hawaii2017/ awsteiner@lit:Hawaii2017
	rsync -Cavzu Hawaii2017Edit/ awsteiner@lit:Hawaii2017Edit

qget:
	cd Hawaii2017; rsync -Cavzun awsteiner@lit:Hawaii2017/* .
	cd Hawaii2017Edit; rsync -Cavzun awsteiner@lit:Hawaii2017Edit/* .

get:
	cd Hawaii2017; rsync -Cavzu awsteiner@lit:Hawaii2017/* .
	cd Hawaii2017Edit; rsync -Cavzu awsteiner@lit:Hawaii2017Edit/* .

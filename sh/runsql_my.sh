#!/bin/sh
#
cwd=`pwd`
#
DBHOST=""
DBPORT="3306"
DBNAME=""
DBUSR=""
DBPW=""
#
MYSQL="mysql"
#
if [ -e "$cwd/.runsql_my" ]; then
	. $cwd/.runsql_my
fi
#
help() {
	echo "$1"
	echo "syntax: `basename $0` [options]"
	echo ""
	echo "  required:"
	echo "        -f FILE ........ SQL file"
	echo "  or"
	echo "        -q QUERY ....... SQL"
	echo "  or"
	echo "        -t ............. show tables"
	echo ""
	echo "  parameters:"
	echo "        -n DBNAME ...... db name [$DBNAME]"
	echo "        -h DBHOST ...... db host [$DBHOST]"
	echo "        -z DBPORT ...... db port [$DBPORT]"
	echo "        -u DBUSR ....... db user [$DBUSR]"
	echo "        -p PW .......... db password"
	echo "  options:"
	echo "        -c ............. TSV output"
	echo "        -v ............. verbose"
	echo ""
	echo "$EXE version: `$MYSQL -V`"
	exit 1
}
#
if [ $# -eq 0 ]; then
	if [ ! -e "$cwd/.runsql_my" ]; then
		echo "Not found: .runsql_my"
		echo "Example:"
		echo "#========cut-here========"
		echo "DBHOST=\"localhost\""
		echo "DBNAME=\"beerdb\""
		echo "#========cut-here========"
	fi
	help ""
fi
#
CSV=""
TSV=""
VERBOSE=""
### Parse options
while getopts f:q:n:h:z:u:p:o:ctv opt ; do
	case "$opt"
	in
	f)      SQLFILE=$OPTARG ;;
	q)      SQL=$OPTARG ;;
	n)      DBNAME=$OPTARG ;;
	h)      DBHOST=$OPTARG ;;
	z)      DBPORT=$OPTARG ;;
	u)      DBUSR=$OPTARG ;;
	p)      DBPW=$OPTARG ;;
	c)      TSV="TRUE" ;;
	t)      SQL="SHOW TABLES" ;;
	v)      VERBOSE="TRUE" ;;
	\?)     help
		exit 1 ;;
	esac
done
#
if [ ! "$DBHOST" -o ! "$DBNAME" ]; then
	help "ERROR: DB specification required."
fi
#
if [ ! "$SQL" -a ! "$SQLFILE" ]; then
	echo "Sql input via -f or -q required."
	help
fi
#
DBOPTS="--no-defaults"
#
if [ "$TSV" ]; then
	DBOPTS="$DBOPTS -ABr"
elif [ "$CSV" ]; then
	DBOPTS="$DBOPTS -ABr"
else
	DBOPTS="$DBOPTS -At"
fi
#
DBOPTS="$DBOPTS --host=$DBHOST --port=$DBPORT --user=$DBUSR"
if [ "$DBPW" ]; then
	DBOPTS="$DBOPTS --password=$DBPW"
fi
#
if [ "$SQLFILE" ]; then
	$MYSQL $DBOPTS $DBNAME <$SQLFILE
else
	$MYSQL $DBOPTS --execute="$SQL" $DBNAME
fi
#

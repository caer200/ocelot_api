# using sphinx

version=0.2
author='Ai, Qianxiang'
projectname='ocelot'

whereisit="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $whereisit
cd ../
sphinx-apidoc . --full -o ./doc -H "$projectname" -A "$author" -V $version
cp $whereisit/conf.py ./doc/
cp $whereisit/index.rst ./doc/
cp quick_start/*.md ./doc/
cp quick_start/*.png ./doc/
cd doc
make html
rm -rf ../docs
mv _build/html ../docs
cd ../
rm -rf doc
touch docs/.nojekyll
cd $whereisit



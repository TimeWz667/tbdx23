
cd dist

git init
git add -A
git commit -m 'deploy'

git push -f git@github.com:TimeWz667/tbdx23.git master:gh-pages

cd ../

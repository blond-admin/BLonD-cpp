#How to clone with features
1. Clone repository
2. Create AppVeyor, TravisCI and Coveralls accounts using your github account
3. Turn on support for BLonD project in each of them
4. Now you have CI for Windows and Linux and Coveralls installed

##Doxygen
1. Create a gh-pages brunch.
2. Create a [Personal Access Token](https://github.com/settings/tokens)
Only enable `public_repo` access for public repositories, "repo" for private.
3. Save the token somewhere as you can only see it once. ([reference](http://stackoverflow.com/a/33125422/1973207))
4. Install the travis gem ([see reference for details](http://stackoverflow.com/a/8113213/1973207)):
```
sudo apt-get install ruby-dev
gem install travis
```
5. Then `cd` into your repository 
6. remove `env->global->secure` in `BLonD-minimal-cpp/.travis.yml`
7. call:
```
travis encrypt GITHUB_API_KEY=<api-token> --add
```
8. update **repository** and **user** names in `git push` at the end of CI script.

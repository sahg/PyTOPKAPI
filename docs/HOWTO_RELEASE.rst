================================
Notes on producing a new release
================================

* Create a release branch off the develop branch, following the model
  proposed in http://nvie.com/posts/a-successful-git-branching-model/

* Generate appropriate source and Windows binary distributions and
  upload to the Github download page.

* Be sure to register the release on PyPI using `python setup.py
  register` - requires authentication.

* Announce release on mailing list.

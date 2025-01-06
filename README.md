# Detection of DNA base modification using nanopore sequencing

In this comprehensive two-day course, participants will learn how to use nanopore sequencing to detect DNA base modifications. Designed to be used to accompany the course and also for self-guided students.
Click [here](https://sathish-t.github.io/nanopore-mod-course) to see the site generated from this github repository.

This repository was forked from a Data Carpentry [repository](https://github.com/datacarpentry/semester-biology)
and is under active development. Some of the text in this README is from the README of that repository.

For the 2024 version of this course website, please go [here](https://tgac.github.io/nanopore-mod-course-2024).
For the 2024 repository, please go [here](https://github.com/tgac/nanopore-mod-course-2024).

## Where is everything

Core teaching materials are stored in `exercises/`, `lectures/`, and 
`materials/`.

Most of the other folders and files support creating the course website using
[Jekyll](http://jekyllrb.com/).

## Creating new pages

If you want to add new exercises, lecture notes, etc. you do this by creating a
[markdown file](http://daringfireball.net/projects/markdown/basics) in the
appropriate directory. Each markdown file needs to start with some information
that tells Jekyll what the page is. This is done using something called YAML,
and the standard YAML for a new exercise would look like this:

```
---
layout: exercise
topic: Topic group of exercise
element: exercise
title: Name of exercise
---
```

This is placed at the very beginning of the markdown file and provides
information on what kind of content it is (e.g., exercise, page, etc.),
the title of the page, and what language it applies to.

The page should then be available at a url based on where the file is located
and what the file name is. So if you created a new exercise in the `exercises/`
folder called `my_awesome_exercise.md` it would be located at:

After pushing to GitHub:
`https://sathish-t.github.io/nanopore-mod-course/exercises/my_awesome_exercise`

## Acknowledgements

Repository forked from [here](https://github.com/datacarpentry/semester-biology).

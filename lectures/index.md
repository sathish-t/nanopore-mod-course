---
layout: page
title: Lectures
---

{% for lecture in site.pages %}
  {% if lecture.layout == 'page' and lecture.element == 'lecture' %}
    * [{{ lecture.title }}]({{ lecture.url | prepend: site.baseurl }})
  {% endif %}
{% endfor %}
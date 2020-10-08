# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Utility functions to create threads running processing functions.
"""

import threading


class Output:
    """ Base class that have an output Queue and basic
    methods to get output items.
    """
    def __init__(self):
        self._done = False
        self._condition = threading.Condition()
        self._tasks = []

    def getTask(self):
        """ This function should be called from a consumer of this
        output instance.
        """
        self._condition.acquire()

        doWait = True
        task = None

        while doWait:
            doWait = False
            if self._tasks:
                task = self._tasks.pop(0)
            elif not self._done:
                self._condition.wait()
                doWait = True

        self._condition.release()

        # Return the task, either None if nothing else should be
        # done, or a task to be processed
        return task

    def _putTask(self, task):
        """ This function should be used by subclasses of Output
        that produces items that will be used by consumers.
        """
        self._condition.acquire()
        self._tasks.append(task)
        self._condition.notify()
        self._condition.release()

    def _setDone(self):
        """ This function should be used by subclasses of Output
        when the output generation is finished. After this getTask
        should return None
        """
        self._condition.acquire()
        self._done = True
        self._condition.notifyAll()
        self._condition.release()


class Generator(Output, threading.Thread):
    def __init__(self, generator):
        threading.Thread.__init__(self)
        Output.__init__(self)
        self._generator = generator

    def run(self):
        for task in self._generator():
            self._putTask(task)

        self._setDone()


class Processor(Generator):
    def __init__(self, processor, inputSource):
        Generator.__init__(self, self._generate)
        self._processor = processor
        self._inputSource = inputSource

    def _generate(self):
        task = self._inputSource.getTask()
        while task is not None:
            yield self._processor(task)  # yield new processed task
            task = self._inputSource.getTask()


class Pipe:
    def __init__(self, generator, *processors):
        self._threads = [Generator(generator)]
        for i, p in enumerate(processors):
            self._threads.append(Processor(p, self._threads[i]))

    def start(self):
        for th in self._threads:
            th.start()

    def join(self):
        for th in self._threads:
            th.join()

    def getTask(self):
        return self._threads[-1].getTask()

